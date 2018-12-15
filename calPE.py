#!/usr/bin/env python3
import sys, os
from glob import glob
import argparse, logging
from datetime import datetime
import pyiast # requires python 3
import numpy as np
import pandas as pd

def wcvf(y, T, P):
    """ Returns the moles of CO2 per kg of adsorbent in the void fraction (wcvf)
    considering the ideal gas law

    :param y: CO2 fraction in the gas phase (-)
    :param T: temperature of the system (K)
    :param P: total pressure (Pa)
    :global vf: void fraction of the bed
    :global ms: mass of adsorbent in per volume of bed (kg/m3)
    """
    R = 8.3144621 # J/mol/K
    return P * y * vf / (R * T * ms)

def ConstructInterpolatorIsothermAtTnew(df, T_0, T_new,
                                        loading_key=None,
                                        pressure_key=None,
                                        hoa_key=None,
                                        fill_value=None):
    """
    Returns InterpolatorIsotherm from pure-component isotherm data measured at
    temperature T0 extrapolated to temperature Tnew using the heat of adsorption.

    :param df: DataFrame Pandas DataFrame with pure-component isotherm tabular
               data measured at temperature T0
    :param T_0: Float temperature at which isotherm in df was measured (T0)
    :param T_new: Float temperature at which we desire to extrapolate the isotherm
               in df
    :param loading_key: String key for loading column in df
    :param pressure_key: String key for pressure column in df
    :param hoa_key: String key for heat of adsorption column in df
    :param fill_value: Float value of loading to assume when an attempt is made
                       to interpolate at a pressure greater than the largest
                       pressure observed in the data

    HOA needs to be in units kJ/mol.
    """

    R = 8.3144621 / 1000.0 # kJ/mol/K
    n = df.shape[0] # number of values in the isotherm
    df_new = pd.DataFrame()
    df_new[pressure_key] = np.zeros((n,))
    df_new[loading_key] = df[loading_key].values

    # for every point, shift pressures according to Classius-Clapyeron eqn
    for i in range(n):
        p_0 = df[pressure_key].iloc[i] # Pa
        HoA = df[hoa_key].iloc[i] # kJ/mol (negative)
        p_new = p_0 * np.exp( HoA / R * (1.0 / T_new - 1.0 / T_0))
        df_new[pressure_key].iloc[i] = p_new

    # return a pyiast.isotherm object
    return pyiast.InterpolatorIsotherm(df_new,
                                       loading_key=loading_key,
                                       pressure_key=pressure_key,
                                       fill_value=fill_value), df_new

def totalQ(struc, Td, Pd):
  """ Returns [0] Q_{thermal} computed from eq. 1 in 10.1039/c4ee02636e (kJ/mol)
              [1] pur, final CO2 purity (-)
              [2] m1, mass CO2 per m3 of bed (kg/m3)
              [3] wct, working capacity CO2 ()
              [4] Qs, sensible heat (kJ/mol)
              [5] Qd, enthalpy of desorption (kJ/mol, defined as positive)
  """
  global ya, yd, qa, wct
  CO2_KGMOL = 0.044 # kg/mol
  # extract or compute the adsorption conditions if not already stored
  if not qa:
    wcvf_a['CO_2'] = wcvf(ya, Ta, Pa)
    wcvf_a['N_2'] = wcvf(1-ya, Ta, Pa)
    # gas uptake @ adsorption
    qa['CO_2'], qa['N_2'] = pyiast.iast(
                                      np.array([ya, 1-ya]) * Pa,
                                      [iso_Ta['CO_2'], iso_Ta['N_2']],
                                      verboseflag=False,
                                      warningoff=True, # no complain about extrapolation
                                     )
    logging.debug("qa['CO_2']: {:.3e}, qa['N_2']: {:.3e}".format(qa['CO_2'], qa['N_2']))
  # Compute desorption conditions
  wcvf_d = {}
  wcvf_d['CO_2'] = wcvf(yd, Td, Pd)
  wcvf_d['N_2'] = wcvf(1-yd, Td, Pd)
  qd, qds = {}, {}
  # gas uptake @ desorption
  qd['CO_2'], qd['N_2'] = pyiast.iast(
                                    np.array([yd, 1-yd]) * Pd,
                                    [iso_Td['CO_2'], iso_Td['N_2']],
                                    verboseflag=False,
                                    warningoff=True, # no complain about extrapolation
                                    adsorbed_mole_fraction_guess=[0.99999999, 0.00000001]
                                   )
  logging.debug("qd['CO_2']: {:.3e}, qd['N_2']: {:.3e}".format(qd['CO_2'], qd['N_2']))
  wc, wct = {}, {}
  for m in ['CO_2', 'N_2']:
    wc[m] = qa[m] - qd[m] # working capacity in the adsorbent (mol_CO2/kg_ADS)
    wct[m] = wc[m] + wcvf_a[m] - wcvf_d[m] # working capacity total = adsorbent + void
  logging.debug("wct['CO_2'] = {:.3e} => {:.3e}".format(wct['CO_2'], wct['CO_2']*CO2_KGMOL))
  if wct['CO_2'] < 0:
    logging.debug("NEGATIVE wct['CO_2']")
    return np.nan, np.nan, np.nan, wct['CO_2'], np.nan, np.nan
  else:
    # Compute CO2 mass, and final purity
    m1 = wct['CO_2'] * ms * CO2_KGMOL
    pur = wct['CO_2'] / (wct['CO_2'] + wct['N_2'])
    logging.debug("m1 = {:.3e}, pur = {:.3e}".format(m1, pur))
    # Compute Q_{thermal} (Qt) = sensible_heat (Qs) + desorption_enthalpy (Qd)
    Qs = cp * ms * (Td - Ta)
    Qd = 0
    for m in ['CO_2', 'N_2']:
      col = iso_df[m].columns[0]
      row = iso_df[m].ix[(iso_df[m][col]-Pd).abs().argsort()[0]]
      Qd += row[iso_df[m].columns[-1]] * wc[m]
    Qd *= ms * 1e3
    Qt = Qs - Qd
    logging.debug("Qs = {:.3e}, Qd = {:.3e}, Qt = {:.3e}".format(Qs, Qd, Qt))
    return Qt, pur, m1, wct['CO_2'], Qs, -Qd

def totalW(Pd, pur):
    """ Compression work (10.1039/c4ee02636e, 'functional representation') """
    plog = -5.4433e5/np.log(1e6/1e3)
    Wp = 7.1723e5+np.log(Pd/1000)*plog
    mp = (2.313-2.102)/(np.sqrt(1e6)-np.sqrt(1e3))
    mpr = 1.102+mp*(np.sqrt(Pd)-np.sqrt(1e3))
    Wt = Wp*(1+mpr*(1./pur-1))
    return Wt

def totalPE(struc, Td, Pd):
  """ Total PE, computed from eq. 2 in 10.1039/c4ee02636e (method = 'carnot') """
  Qt, pur, m1, wct, Qs, Qd = totalQ(struc, Td, Pd) # Compute mainly the PE and CO2 purity
  if wct < 0: # Negative working capacity
      return np.nan, np.nan, np.nan, np.nan, wct, np.nan, np.nan, np.nan, np.nan
  Wt = totalW(Pd, pur) # Compute compression work
  logging.debug("Wt = {:.3e}".format(Wt))
  if method == "carnot":
      nu = (Td + 10 - 283)/(Td + 10)
      logging.debug("nu = {:.3e}".format(nu))
      Qteff = 0.75 * Qt * nu / m1
      Qs =    0.75 * Qs * nu / m1
      Qd =    0.75 * Qd * nu / m1
  elif method == "linear":
      Qteff = Qt * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
      Qs = Qs * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
      Qd = Qd * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
  logging.debug("Qt = {:.3e}, Qteff = {:.3e}".format(Qt, Qteff))
  PE = Qteff + Wt
  return PE, Qteff, Wt, m1, wct, Qs, Qd, pur, Qt

def main(args):
    global cp, ms, ya, yd, Ta, Pa, vf, gas, totE, pCO2, method
    global iso_Ta, iso_df, T_df, iso_Td, qa, wcvf_a
    # Read parameters
    method = args.method
    vf = args.vf
    cp = args.cp
    yd = args.yd
    ms = args.rho * ( 1. - vf ) # Mass (kg) of adsorbent in per m3 of bed
    logging.debug("ms = {:.3e}".format(ms))
    if args.comp == 'coal':
        totE = 6631.2
        pCO2 = 1.80 # kg_CO2/kg_coal (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.14 # 14:86 ratio
        Ta = 313.0
        Pa = 101325.0
    elif args.comp == 'ng':
        totE = 21023.26
        pCO2 = 3.23 # kg_CO2/kg_ng (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.04 # 4% CO2
        Ta = 313.0
        Pa = 101325.0
    elif args.comp == 'air':
        ya = 0.0004 # 400ppm CO2
        Ta = 288.0
        Pa = 101325.0
    T_df = {}
    iso_df = {}
    iso_Ta = {}
    qa = {}
    wcvf_a = {}
    for m in ['CO_2','N_2']:
        df_dir = os.path.join(args.datapath, args.struc, m) # dir containing isotherms
        df_file = glob(os.path.join(df_dir, '*.csv'))[0] # .csv files for the isotherms
        T_df[m] = int(os.path.splitext(os.path.basename(df_file))[0][:-1]) # temperature of the isotherms
        iso_df[m] = pd.read_csv(df_file, sep=' ') #values of the isotherm
        # Extrapolate the pure-component isotherm at Ta
        iso_Ta[m], _ = ConstructInterpolatorIsothermAtTnew(
                              iso_df[m], T_df[m], Ta,
                              loading_key="loading(mol/kg)",
                              pressure_key="pressure(Pa)",
                              hoa_key="HoA(kJ/mol)",
                              fill_value=iso_df[m]["loading(mol/kg)"].max())
    # Set the range and the step for the Td and Pd, tested to find the min PE:
    # the range was changed from 10.1039/c4ee02636e, pag. 4136, 2nd col.
    Td_min = Ta + 20.0
    Td_step = 10.0
    Td_max = Td_min + 100.0
    Pd_min = 0.01 * 101325.0
    Pd_step = 0.02 * 101325.0
    Pd_max = 1.5 * 101325.0
    if args.process == 'PSA':
        Td_range = np.array([Td_min])
        Pd_range = np.arange(Pd_min,Pd_max,Pd_step)
    elif args.process == 'TSA':
        Td_range = np.arange(Td_min,Td_max,Td_step)
        Pd_range = np.array([Pa])
    elif args.process == 'TPSA':
        Td_range = np.arange(Td_min,Td_max,Td_step)
        Pd_range = np.arange(Pd_min,Pd_max,Pd_step)
    # Compute the PE at different Td, Pd
    data = [] # collect all the data at different Pd, Td
    iso_Td = {}
    for Td in Td_range:
        for m in ['CO_2','N_2']:
            # Extrapolate the pure-component isotherm at Td
            iso_Td[m], _ = ConstructInterpolatorIsothermAtTnew(
                                 iso_df[m], T_df[m], Td,
                                 loading_key="loading(mol/kg)",
                                 pressure_key="pressure(Pa)",
                                 hoa_key="HoA(kJ/mol)",
                                 fill_value=iso_df[m]["loading(mol/kg)"].max())
        for Pd in Pd_range:
            # Compute the PE @ Td and Pd
            logging.debug("**** Evaluating: Td={:.1f}, Pd={:.3f} *****".format(Td, Pd))
            PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt = totalPE(args.struc, Td, Pd)
            if wct > 0 : #Positive working capacity
                logging.debug("PE = {:.3e}".format(PE))
                data.append([Td, Pd, PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt])
    # Find the conditions for the min PE and extract the results
    if len(data) == 0: # All wc are negative
        return "{:s}: Unfeasible process!".format(args.struc)
    data = np.array(data)
    data_minPE = data[np.argmin(data.T[2])]
    logging.debug("data_minPE:")
    logging.debug(data_minPE)
    finPE = data_minPE[2]/1e6   # minimum PE (MJ/kg)
    finP = data_minPE[1]/101325.0 # desorption pressure (bar)
    finT = data_minPE[0]          # desorption temperature (K)
    if args.comp in ['coal','ng']:
        finEL = finPE*1e3 * pCO2 / totE # fraction of electricity loss (-)
    elif args.comp == 'air':
        finEL = np.nan
    finCap = data_minPE[3]/1e3  # heat requirement (kJ/kgCO2)
    finComp = data_minPE[4]/1e3 # compression work (kJ/kgCO2)
    finMProd = data_minPE[5]    # Mass of CO2 produced (kg)
    finWC = data_minPE[6]       # CO2 working capacity (mol/kg)
    finPur = data_minPE[9]      # fraction of CO2 purity (-)
    finQs = data_minPE[7]/1e3   # (not printed) Heat required to carry out the separation, Q_thermal (kJ/kgCO2)
    finQd = data_minPE[8]/1e3   # (not printed) Energy to heat the adsorbent for desorption, sensible heat (kJ/kgCO2)
    finQt = data_minPE[10]/1e3  # (not printed) Energy to undo the adsorption process, Dh (kJ/kgCO2)
    return "{:s}: PE(MJ/kg)= {:.3f} Pd(bar)= {:.3f} Td(K)= {:.1f} ".format(
            args.struc,finPE,finP,finT) + \
           "EL(-) = {:.3f} Q(kJ/kg)= {:.2f} Wcomp(kJ/kg)= {:.2f} ".format(
            finCap,finComp,finEL) + \
           "M(kg)= {:.3f} WC(mol/kg)= {:.3f} pur(-)= {:.3f}".format(
            finMProd,finWC,finPur)


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Program to compute the parasitic energy from single component isotherms",
    formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument(
          "struc",
          help="Name of the adsorbent framework.")
  parser.add_argument(
          "rho",
          type=float,
          help="Density of the adsorbent framework (kg/m3)")
  parser.add_argument(
          "comp",
          choices=["coal", "ng", "air"],
          help="Compositon of the mixture with CO2.")
  parser.add_argument(
          "-vf",
          type=float,
          dest="vf",
          default=0.35,
          help="Void fraction of the adsorbent bed (default: 0.35).")
  parser.add_argument(
          "-process",
          dest="process",
          default="TPSA",
          choices=["TPSA", "TSA", "PSA"],
          help="Process used for the CO2 sequestration.")
  parser.add_argument(
          "-cp",
          type=float,
          dest="cp",
          default=985.0,
          help="Choice for the Cp of the adsorbent:\n" +
               "for nanoporous materials it should range between 761.0\n" +
               "and 1210.0 J/kg/K (default: 985.0 J/Kg/K, the average)")
  parser.add_argument(
          "-yd",
          type=float,
          dest="yd",
          default=0.99,
          help="Required output CO2 fraction at desorption (default: 0.99)")
  parser.add_argument(
          "-method",
          dest="method",
          default="carnot",
          choices=["carnot", "linear"],
          help="Method to compute the electricity conversion efficiency.")
  parser.add_argument(
          "-datapath",
          type=str,
          dest="datapath",
          default="./test/",
          help="Path containing the isotherms in .csv format (default: ./ccsdata)")
  parser.add_argument(
         "-l",
         "--log",
         action = "store_true",
         help="Write .log file for debug")
  args = parser.parse_args()
  if args.log:
    # Set up the debug .log file
    now = datetime.now().isoformat()
    logfile = 'calPE_%s_%s.log' % (args.struc, now)
    logging.basicConfig(
      filename=logfile, format='%(message)s', level=logging.DEBUG
    )
  # Write input in .log, run the main program and print results in .log and on screen
  logging.debug(args)
  results = main(args)
  logging.debug(results)
  print(results)
