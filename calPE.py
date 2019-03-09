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
              [1] Qs, sensible heat (kJ/mol)
              [2] Qd, enthalpy of desorption (kJ/mol, defined as positive)
              [3] WCv, mass CO2 per m3 of bed (kg/m3)
              [4] pur, final CO2 purity (-)

  """
  global ya, yd, qa, wct
  CO2_KGMOL = 0.044 # kg/mol
  # extract or compute the adsorption conditions if not already stored
  if not qa:
    logging.debug("Ta = {:.3e}, Pa = {:.3e}, ya = {:.3e}, yd = {:.3e}".format(Ta,Pa,ya,yd))
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
  WCv = wct['CO_2'] * ms * CO2_KGMOL
  logging.debug("wct['CO_2'] = {:.3e} => WCv = {:.3e}".format(wct['CO_2'], WCv))
  if WCv < 0:
    logging.debug("NEGATIVE wct['CO_2']")
    return np.nan, np.nan, np.nan, WCv, np.nan
  else:
    # Compute volumetric working capacity, and final purity

    pur = wct['CO_2'] / (wct['CO_2'] + wct['N_2'])
    logging.debug("WCv = {:.3e}, pur = {:.3e}".format(WCv, pur))
    # Compute Q_{thermal} (Qt) = sensible_heat (Qs) + desorption_enthalpy (Qd)
    Qs = cp * ms * (Td - Ta)
    Qd = 0
    for m in ['CO_2', 'N_2']:
      col = iso_df[m].columns[0]
      row = iso_df[m].ix[(iso_df[m][col]-Pd).abs().argsort()[0]]
      Qd += -(row[iso_df[m].columns[-1]] * wc[m])
    Qd *= ms * 1e3
    Qt = Qs + Qd # Note that Qd is negative!
    logging.debug("Qs = {:.3e}, Qd = {:.3e}, Qt = {:.3e}".format(Qs, Qd, Qt))
    return Qt, Qs, Qd, WCv, pur

def totalWcomp(Pd, pur):
    """ Compression work (kJ/kgCO2) from P_desorption (Pa) to 150bar. This uses
    a 'functional representation' as described in 10.1039/c4ee02636e.
    """
    plog = -5.4433e5/np.log(1e6/1e3)
    Wp = 7.1723e5+np.log(Pd/1000)*plog
    mp = (2.313-2.102)/(np.sqrt(1e6)-np.sqrt(1e3))
    mpr = 1.102+mp*(np.sqrt(Pd)-np.sqrt(1e3))
    Wcomp = Wp*(1+mpr*(1./pur-1))
    return Wcomp

def totalPE(struc, Td, Pd):
  """ Total PE, computed from eq. 2 in 10.1039/c4ee02636e (eleff = 'carnot') """
  Qt, Qs, Qd, WCv, pur = totalQ(struc, Td, Pd) # Compute mainly the PE and CO2 purity
  if WCv < 0: # Negative working capacity: return NaN for all values but WCv
      return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, WCv, np.nan
  # Compute compression work from Pd to 150bar
  Wcomp = totalWcomp(Pd, pur)
  logging.debug("Wcomp = {:.3e}".format(Wcomp))
  if eleff == "carnot":
      nu = (Td + 10 - 283)/(Td + 10)
      logging.debug("nu = {:.3e}".format(nu))
      Qteff = 0.75 * Qt * nu / WCv
      Qseff = 0.75 * Qs * nu / WCv
      Qdeff = 0.75 * Qd * nu / WCv
  elif eleff == "linear":
      Qteff = Qt * (-0.08037 + 0.002326 * (Td - 273 + 10)) / WCv
      Qseff = Qs * (-0.08037 + 0.002326 * (Td - 273 + 10)) / WCv
      Qdeff = Qd * (-0.08037 + 0.002326 * (Td - 273 + 10)) / WCv
  logging.debug("Qt = {:.3e}, Qteff = {:.3e}".format(Qt, Qteff))
  PE = Qteff + Wcomp
  return PE, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur

def mainPE(struc, gasin, rho, vf, process, cp, yd, eleff, opt, T_iso, iso_df):
    """Main script to compute CO2 parasitic energy from single component
    CO2 and N2 isotherms

    :str struc: name of the adsorbent
    :str gasin: choice of the input mixture
    :float rho: density of the adsorbent (kg/m3)
    :float vf:  void fraction of the bed, e.g., 0.35
    :str process: choice of the process, e.g., TPSA
    :float cp: heat capacity of the adsorbent (J/kg/K), e.g., 985.0
    :float yd: target CO2 purity, e.g., 0.99
    :str eleff: choice of the method to compute electr. eff., e.g., 'carnot'
    :str opt: choice of the parameter to optimize for the process, e.g., 'PE'
    :[float,float] T_iso: temperature of the used isotherms for CO2 and N2
    :[df,df] iso_df: dataframe with isotherms for CO2 and N2, see examples
    """
    global cp, ms, ya, yd, Ta, Pa, vf, gas, totE, eleff
    global iso_Ta, iso_df, T_df, iso_Td, qa, wcvf_a

    # Initialize calculation
    ms = rho * ( 1. - vf ) # Mass (kg) of adsorbent in per m3 of bed

    if gasin == 'coal':
        totE = 6631.2 / 1.80 # kJ/kg_CO2 (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.14 # 14:86 molar ratio
        Ta = 313.0
        Pa = 101325.0
    elif gasin == 'ng':
        totE = 21023.26 / 3.23  # kJ/kg_CO2 (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.04 # 4% molar CO2
        Ta = 313.0
        Pa = 101325.0
    elif gasin == 'air':
        totE = np.nan # ill defined
        pCO2 = np.nan # ill defined
        ya = 0.0004 # 400ppm molar CO2
        Ta = 288.0
        Pa = 101325.0
    iso_Ta = {}
    for m in ['CO_2','N_2']:
        # Extrapolate the pure-component isotherm at Ta
        iso_Ta[m], _ = ConstructInterpolatorIsothermAtTnew(
                              iso_df[m], T_df[m], Ta,
                              loading_key="loading(mol/kg)",
                              pressure_key="pressure(Pa)",
                              hoa_key="HoA(kJ/mol)",
                              fill_value=iso_df[m]["loading(mol/kg)"].max())

    # Set the range and the step for the Td and Pd, tested to find the min PE:
    # the range was changed from {10.1039/c4ee02636e, pag. 4136, 2nd col},
    # to have less points and be faster.
    Td_min = Ta + 20.0
    Td_step = 10.0
    Td_max = Td_min + 100.0
    Pd_min = 0.01 * 101325.0
    Pd_step = 0.02 * 101325.0
    Pd_max = 1.5 * 101325.0
    if process == 'PSA':
        Td_range = np.array([Td_min])
        Pd_range = np.arange(Pd_min,Pd_max,Pd_step)
    elif process == 'TSA':
        Td_range = np.arange(Td_min,Td_max,Td_step)
        Pd_range = np.array([Pa])
    elif process == 'TPSA':
        Td_range = np.arange(Td_min,Td_max,Td_step)
        Pd_range = np.arange(Pd_min,Pd_max,Pd_step)

    # Compute the PE at different Td, Pd
    data = [] # collect all the data at different Pd, Td
    qa = {}
    wcvf_a = {}
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
            PE, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur = totalPE(args.struc, Td, Pd)
            if WCv > 0 : #Positive working capacity
                logging.debug("PE = {:.3e}".format(PE))
                data.append([PE, Pd, Td, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur])

    # Find the optimal conditions required and combine the results in a dict
    if len(data) == 0: # All wc are negative
        results_dict = {'process_feasible' : False}
        return results_dict
    data = np.array(data)
    if opt == "PE":
        data_opt = data[np.argmin(data.T[0])]
    elif opt == "Q":
        data_opt = data[np.argmin(data.T[3])]
    elif opt == "WC":
        data_opt = data[np.argmax(data.T[8])]
    elif opt == "pur":
        data_opt = data[np.argmax(data.T[9])]
    logging.debug("data_opt:")
    logging.debug(data_opt)
    results_dict = {}
    results_dict['process_feasible'] = True
    results_dict['PE'] = data_opt[0]/1e6
    results_dict['PE_units'] = 'MJ/kg'
    results_dict['PE_descr'] = 'parassitic energy at optimal process conditions'
    results_dict['P'] = data_opt[1]/101325.0
    results_dict['P_units'] = 'bar'
    results_dict['PE_descr'] = 'desorption pressure'
    results_dict['T'] = data_opt[2]
    results_dict['T_units'] = 'K'
    results_dict['T_descr'] = 'desorption temperature'
    results_dict['eloss'] = finPE*1e3 / totE
    results_dict['eloss_units'] = 'kJ/kJ'
    results_dict['eloss_descr'] = 'fraction of electricity loss (np.nan for "air")'
    results_dict['Qteff'] = data_opt[3]/1e6
    results_dict['Qteff_units'] = 'MJ/kgCO2'
    results_dict['Qteff_descr'] = 'heat requirement'
    results_dict['Qt'] = data_opt[6]/1e6
    results_dict['Qt_units'] = 'MJ/kgCO2'
    results_dict['Qt_descr'] = 'Not-effective heat requirement'
    results_dict['Qs'] = data_opt[4]/1e6
    results_dict['Qs_units'] = 'MJ/kgCO2'
    results_dict['Qs_descr'] = 'Sensible heat'
    results_dict['Qd'] = data_opt[5]/1e6
    results_dict['Qd_units'] = 'MJ/kgCO2'
    results_dict['Qd_descr'] = 'desorption enthalpy'
    results_dict['Wcomp'] = data_opt[7]/1e6
    results_dict['Wcomp_units'] = 'MJ/kgCO2'
    results_dict['Wcomp_descr'] = 'compression work'
    results_dict['WCv'] = data_opt[8]
    results_dict['WCv_units'] = 'kgCO2/m3'
    results_dict['WCv_descr'] = 'volumetric working capacity'
    results_dict['WCg'] = data_opt[8]/ms
    results_dict['WCg_units'] = 'kgCO2/kg'
    results_dict['WCg_descr'] = 'gravimetric working capacity'
    results_dict['Pur'] = data_opt[9]
    results_dict['Pur_units'] = 'mol/mol
    results_dict['Pur_descr'] = 'fraction of CO2 purity'
    return results_dict

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Program to compute the parasitic energy from single component isotherms",
    formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument(
          "struc",
          help="Name of the adsorbent framework.")
  parser.add_argument(
          "gasin",
          choices=["coal", "ng", "air"],
          help="Compositon of the input gas mixture, i.e. post-combustion\n" +
               "flue gas from coal or natural gas, or air.")
  parser.add_argument(
          "-rho",
          type=float,
          dest="rho",
          default=None,
          help="Density of the adsorbent framework (kg/m3).\n" +
               "(default: value readen from datapath/struc/rho.csv)")
  parser.add_argument(
          "-vf",
          type=float,
          dest="vf",
          default=0.35,
          help="Void fraction of the adsorbent bed.\n" +
               "(default: 0.35)")
  parser.add_argument(
          "-process",
          dest="process",
          default="TPSA",
          choices=["TPSA", "TSA", "PSA"],
          help="Process used for the CO2 sequestration.\n" +
               "(default: TPSA)")
  parser.add_argument(
          "-cp",
          type=float,
          dest="cp",
          default=None,
          help="Choice for the heat adsorption of the adsorbent:\n" +
               "for nanoporous materials it should range between 761.0\n" +
               "and 1210.0 J/kg/K.\n" +
               "(default: readen from datapath/struc/cp.csv if present,\n" +
               "          otherwise 985.0 J/kg/K, the average)")
  parser.add_argument(
          "-yd",
          type=float,
          dest="yd",
          default=0.99,
          help="Required output CO2 fraction at desorption.\n" +
               "(default: 0.99)")
  parser.add_argument(
          "-eleff",
          dest="eleff",
          default="carnot",
          choices=["carnot", "linear"],
          help="Method to compute the electricity conversion efficiency.\n" +
               "(default: carnot)")
  parser.add_argument(
          "-opt",
          dest="opt",
          default="PE",
          choices=["PE", "Q", "WC", "pur"],
          help="Optimizes for the optimal condition: lowest PE or Q,\n" +
               "or highest working capacity (WC) or pur(ity).\n" +
               "(default: PE)")
  parser.add_argument(
          "-datapath",
          type=str,
          dest="datapath",
          default="./test/",
          help="Path containing the isotherms in .csv format.\n" +
               "(default: ./test/)")
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
  logging.debug(args)

  # Convert varables into the arguments to feed the mainPE function
  if args.rho == None:
        with open(os.path.join(args.datapath, args.struc, "rho.csv")) as f:
            args.rho = float(f.readline().strip()) #density in kg/m3

  if args.cp == None:
        try:
            with open(os.path.join(args.datapath, args.struc, "cp.csv")) as f:
                cp = float(f.readline().strip()) #cp in J/kg/K
        except:
            cp = 948.0 #average for MOFs
  else:
        cp = args.cp

  for m in ['CO_2','N_2']:
    iso_dir = os.path.join(args.datapath, args.struc, m) # dir containing isotherms
    iso_file = glob(os.path.join(df_dir, '*.csv'))[0] # .csv files for the isotherms
    T_iso[m] = int(os.path.splitext(os.path.basename(df_file))[0][:-1]) # read temperature of the isotherm from the filename
    iso_df[m] = pd.read_csv(df_file, sep=' ') #values of the isotherm

  # Run the main function
  res = mainPE(struc=args.struc,
               gasin=args.gasin,
               rho=args.rho,
               vf=args.vf,
               process=args.process,
               cp=args.cp,
               yd=args.yd,
               eleff=args.eleff,
               opt=args.eleff,
               T_iso=T_iso,
               iso_df=iso_df,
              )

  # Convert results into a single-line string and print it
  results_str="{:s}: ".format(args.struc)
  results_str+="PE(MJ/kg)= {:.3f}: ".format(res['PE'])
  results_str+="Pd(bar)= {:.2f} ".format(res['P'])
  results_str+="Td(K)= {:.1f} ".format(res['T'])
  results_str+="EL(J/J) = {:.3f} ".format(res['EL'])
  results_str+="Q(MJ/kg)= {:.3f} ".format(res['Qteff'])
  results_str+="Wcomp(MJ/kg)= {:.3f} ".format(res['Wcomp'])
  results_str+="WCv(kg/m3)= {:.3f} ".format(res['WCv'])
  results_str+="WCg(kg/kg)= {:.3f} ".format(res['WCg'])
  results_str+="pur(mol/mol)= {:.3f}".format(res['Pur'])
  logging.debug(results_str)
  print(results_str)
