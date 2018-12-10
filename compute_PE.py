#!/usr/bin/env python3

import sys, math, argparse, logging, os
from datetime import datetime
import pyiast # requires python 3
import numpy
import ccsdb, utils
from fnmatch import fnmatch
from glob import glob
import pandas as pd

def wcvf(y, T, P):
  """ Compute the mol of CO2 per kg of adsorbent (wcvf) considering the ideal gas law """
  return P * y * vf / (utils.R * T * ms)

def totalQ(struc, Td, Pd):
  """ Q_{thermal} computed from eq. 1 in 10.1039/c4ee02636e """
  global ya, yd, qa, wct
  # extract or compute the adsorption conditions if not already stored
  if not qa:
    wcvf_a['CO_2'] = wcvf(ya, Ta, Pa)
    wcvf_a['N_2'] = wcvf(1-ya, Ta, Pa)
    if model == 'IAST':
      qa['CO_2'], qa['N_2'] = pyiast.iast(
                                      numpy.array([ya, 1-ya]) * Pa,
                                      [iso_Ta['CO_2'],
                                      iso_Ta['N_2']],
                                      verboseflag=False,
                                     )
      logging.debug("qa['CO_2']: %r, qa['N_2']: %r" % (qa['CO_2'], qa['N_2']))
    elif model == 'COMP':
        kqv = {'CO_2':[1, 0, 0, 1, 1, 0, 0, 1],
               'N_2': [1, 0, 0, 1, 1, 0, 0, 1]}
        qa['CO_2'], qa['N_2'] = utils.compQ(Ta, Pa, kqv, ya, dH)
        logging.debug("qa['CO_2']: %r, qa['N_2']: %r" % (qa['CO_2'], qa['N_2']))
  # Compute desorption conditions
  wcvf_d = {}
  wcvf_d['CO_2'] = wcvf(yd, Td, Pd)
  wcvf_d['N_2'] = wcvf(1-yd, Td, Pd)
  qd, qds = {}, {}
  if model == 'IAST':
    qd['CO_2'], qd['N_2'] = pyiast.iast(
                                    numpy.array([yd, 1-yd]) * Pd,
                                    [iso_Td['CO_2'],
                                    iso_Td['N_2']],
                                    verboseflag=False,
                                    adsorbed_mole_fraction_guess=[0.99999999, 0.00000001]
                                   )
    logging.debug("qd['CO_2']: %r, qd['N_2']: %r" % (qd['CO_2'], qd['N_2']))
  elif model == 'COMP':
      qd['CO_2'], qd['N_2'] = utils.compQ(Td, Pd, kqv, yd, dH)
      logging.debug("qd['CO_2']: %r, qd['N_2']: %r" % (qd['CO_2'], qd['N_2']))
  wc, wct = {}, {}
  for m in ['CO_2', 'N_2']:
    wc[m] = qa[m] - qd[m]
    wct[m] = wc[m] + wcvf_a[m] - wcvf_d[m]
  if wc['CO_2'] < 0:
    return -1, -1, 1, 1, -1, -1
  else:
    Qs = cp * ms * (Td - Ta)
    Qd = 0
    for m in ['CO_2', 'N_2']:
      col = iso_df[m].columns[0]
      row = iso_df[m].ix[(iso_df[m][col]-Pd).abs().argsort()[0]]
      local_dH = row[iso_df[m].columns[-1]]
      Qd += -local_dH * wc[m]
    Qd *= ms * 1e3
    pur = wct['CO_2'] / (wct['CO_2'] + wct['N_2'])
    m1 = wct['CO_2'] * ms * 0.044
    logging.debug("Qs = %r, Qd = %r" % (Qs, Qd))
    logging.debug("wct['CO_2'] = %r => %r" % (wct['CO_2'], wct['CO_2']*0.044))
    logging.debug("m1 = %f, pur = %f" % (m1, pur))
    Qt = Qs - Qd
    return Qt, pur, m1, wct['CO_2'], Qs, -Qd

def totalW(Pd, pur):
    """ Compression work (10.1039/c4ee02636e, 'we developed a functional representation') """
    plog = -5.4433e5/math.log(1e6/1e3)
    Wp = 7.1723e5+math.log(Pd/1000)*plog
    mp = (2.313-2.102)/(math.sqrt(1e6)-math.sqrt(1e3))
    mpr = 1.102+mp*(math.sqrt(Pd)-math.sqrt(1e3))
    Wt = Wp*(1+mpr*(1./pur-1))
    return Wt

def totalPE(struc, Td, Pd):
  """ Total PE, computed from eq. 2 in 10.1039/c4ee02636e (method = 'carnot') """
  Qt, pur, m1, wct, Qs, Qd = totalQ(struc, Td, Pd) # Compute mainly the PE and CO2 purity
  if Qt < 0: return -1, -1, -1, -1, -1, -1, -1, -1, -1
  Wt = totalW(Pd, pur) # Compute compression work
  logging.debug("Wt = %f" % Wt)
  if method == "carnot":
      nu = (Td + 10 - 283)/(Td + 10)
      logging.debug("nu = %r" % nu)
      Qteff = 0.75 * Qt * nu / m1
      Qs =    0.75 * Qs * nu / m1
      Qd =    0.75 * Qd * nu / m1
  elif method == "linear":
      Qteff = Qt * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
      Qs = Qs * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
      Qd = Qd * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
  logging.debug("Qt = %r, Qteff = %r" % (Qt, Qteff))
  PE = Qteff + Wt
  return PE, Qteff, Wt, m1, wct, Qs, Qd, pur, Qt

def main(args):
    global cp, ms, ya, yd, Ta, Pa, dH, vf, gas, model, totE, pCO2, method
    global iso_Ta, iso_df, T_df, iso_Td, qa, wcvf_a
    # Read parameters
    model = args.model
    method = args.method
    vf = args.vf
    cp = args.cp
    yd = args.yd
    ms = args.rho * ( 1. - vf ) # Mass (kg) of adsorbent in per m3 of bed
    logging.debug('ms = %f' %  ms)
    dH = {'CO_2': [-1, -1], 'N_2':[-1, -1]}
    if args.comp == 'coal':
        totE = 6631.2
        pCO2 = 1.80
        ya = 0.14
        Ta = 313.0
        Pa = 101325.0
        Td_min = 333.0
    elif args.comp == 'NG':
        totE = 21023.26
        pCO2 = 3.23
        ya = 0.04
        Ta = 313.0
        Pa = 101325.0
        Td_min = 333.0
    elif args.comp == 'air':
        ya = 0.004
        Ta = 288.0
        Pa = 101325.0
        Td_min = 308.0
    T_df = {}
    iso_df = {}
    iso_Ta = {}
    for m in ['CO_2','N_2']:
        df_dir = os.path.join(args.datapath, args.struc, m) # dir containing isotherms
        df_file = glob(os.path.join(df_dir, '*.csv'))[0] # .csv files for the isotherms
        T_df[m] = int(os.path.splitext(os.path.basename(df_file))[0][:-1]) # temperature of the isotherms
        iso_df[m] = pd.read_csv(df_file, sep=' ') #values of the isotherm
        # Extrapolate the pure-component isotherm at Ta
        iso_Ta[m], _ = utils.ConstructInterpolatorIsothermAtTnew(
                              iso_df[m], T_df[m], Ta,
                              loading_key="loading(mol/kg)",
                              pressure_key="pressure(Pa)",
                              hoa_key="HoA(kJ/mol)",
                              fill_value=iso_df[m]["loading(mol/kg)"].max())
    # Set the range and the step for the Td and Pd, tested to find the min PE:
    # the range is coherent with 10.1039/c4ee02636e, pag. 4136, 2nd col.
    Td_step = 10.0
    Td_max = Td_min+100*Td_step
    Pd_min = 0.01*101325.0
    Pd_step = 0.01*101325.0
    Pd_max = 3.0*101325.0
    if args.process == 'PSA':
        Td_range = numpy.array([Td_min])
        Pd_range = numpy.arange(Pd_min,Pd_max,Pd_step)
    elif args.process == 'TSA':
        Td_range = numpy.arange(Td_min,Td_max,Td_step)
        Pd_range = numpy.array([1.0*101325.0])
    elif args.process == 'TPSA':
        Td_range = numpy.arange(Td_min,Td_max,Td_step)
        Pd_range = numpy.arange(Pd_min,Pd_max,Pd_step)
    # Compute the PE at different Td, Pd
    data = [] # collect all the data at different Pd, Td
    iso_Td = {}
    qa = {}
    wcvf_a = {}
    for Td in Td_range:
        for m in ['CO_2','N_2']:
            # Extrapolate the pure-component isotherm at Td
            iso_Td[m], _ = utils.ConstructInterpolatorIsothermAtTnew(
                                 iso_df[m], T_df[m], Td,
                                 loading_key="loading(mol/kg)",
                                 pressure_key="pressure(Pa)",
                                 hoa_key="HoA(kJ/mol)",
                                 fill_value=iso_df[m]["loading(mol/kg)"].max())
        for Pd in Pd_range:
            # Compute the PE @ Td and Pd
            PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt = totalPE(args.struc, Td, Pd)
            logging.debug("Td=%r, Pd=%r, PE=%r" %(Td, Pd, PE))
            if PE > 0: #Needed because errors give negative PE
                data.append([Td, Pd, PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt])
    # Find the conditions for the min PE and extract the results
    data = numpy.array(data)
    data_minPE = data[numpy.argmin(data.T[2])]
    logging.debug("data_minPE: %r" % data_minPE)
    finPE = round(data_minPE[2]/1000., 6)   # minimum PE (kJ/kg)
    finP = round(data_minPE[1]/101325.0, 3) # desorption pressure (bar)
    finT = round(data_minPE[0], 1)          # desorption temperature (K)
    if args.comp in ['coal','NG']:
        finEL = round(finPE * pCO2 / totE, 5) # fraction of electricity loss (-)
    elif args.comp == 'air':
        finEL = -1
    finCap = round(data_minPE[3]/1000., 6)  # heat requirement (kJ/kgCO2)
    finComp = round(data_minPE[4]/1000., 6) # compression work (kJ/kgCO2)
    finMProd = round(data_minPE[5], 5)      # Mass of CO2 produced (kg)
    finWC = round(data_minPE[6], 5)         # CO2 working capacity (mol/kg)
    finPur = round(data_minPE[9], 5)        # fraction of CO2 purity (-)
    finQs = round(data_minPE[7]/1000., 6)   # (not printed) Heat required to carry out the separation, Q_thermal (kJ/kgCO2)
    finQd = round(data_minPE[8]/1000., 6)   # (not printed) Energy to heat the adsorbent for desorption, sensible heat (kJ/kgCO2)
    finQt = round(data_minPE[10]/1000., 6)  # (not printed) Energy to undo the adsorption process, Dh (kJ/kgCO2)
    # Print the results in the log, on screen and in the database.yml
    finPE_id = [finPE, finP, finT, finEL, finCap, finComp, finMProd, finWC, finPur]
    logging.debug("%s %r" %  (args.struc, finPE_id))
    print (args.struc, finPE_id)
    # END of main

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
          choices=["coal", "NG", "air"],
          help="Compositon of the mixture with CO2.")
  parser.add_argument(
          "-model",
          dest="model",
          default="IAST",
          choices=["IAST", "COMP"],
          help="Model for the mixture adsorption:\n" +
               "IAST (Ideal Adsorbed Solution Theory) (default)\n" +
               "COMP (COMPetititive Langmuir Adsorption model)")
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
          default="./ccsdata/",
          help="Path containing the isotherms in .csv format (default: ./ccsdata)")
  parser.add_argument(
         "-l",
         "--log",
         action = "store_true",
         help="Write .log file for debug")
  args = parser.parse_args()
  if args.log:
    now = datetime.now().isoformat()
    logfile = 'calPE_%s_%s.log' % (args.struc, now)
    logging.basicConfig(
      filename=logfile, format='%(message)s', level=logging.DEBUG
    )
  main(args)
