#!/usr/bin/env python3

import sys, math, argparse, logging, os
from datetime import datetime
import pyiast # requires python 3
import numpy
import ccsdb, utils, deltaT
from fnmatch import fnmatch
from glob import glob
import pandas as pd

cp, ms, Ta, Pa, vf, gas = 0, 0, 0, 0, 0, 0
dH, qa, qas, kqv, wcvf_a = {}, {}, {}, {}, {}
iso_ta, iso_df, T_df, iso_td = {}, {}, {}, {}
model = ''
y1 = []

db = 'coreDB.yml'
struc = 'ABAVIJ'
comp = 'coal'
model = 'IAST'
limit = 'IAST'

def totalPE(db, struc, Td, Pd):
  Qt, pur, m1, wct, Qs, Qd = totalQ(db, struc, Td, Pd)
  if Qt < 0: return -1, -1, -1, -1, -1, -1, -1, -1, -1
  Wt = totalW(Pd, pur)
  logging.debug("Wt = %f" % Wt)
  e = (Td + 10 - 283)/(Td + 10)
  logging.debug("e = %r" % e)
  if method == "carnot":
    Qteff = 0.75 * Qt * e / m1
    Qs = Qs * 0.75 * e / m1
    Qd = Qd * 0.75 * e / m1
  else:
    Qteff = Qt * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
    Qs = Qs * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
    Qd = Qd * (-0.08037 + 0.002326 * (Td - 273 + 10)) / m1
  logging.debug("Qt = %r, Qteff = %r" % (Qt, Qteff))
  return Qteff + Wt, Qteff, Wt, m1, wct, Qs, Qd, pur, Qt

def wcvf(y, T, P):
  return P * y * vf / (utils.R * T * ms)

def totalQ(db, struc, Td, Pd):
  global y1, qa, qas, wct
  if not y1: y1 = ccsdb.getHead(db, 'IastY1')
  # adsorption conditions if not done yet
  if not qa:
    wcvf_a['CO_2'] = wcvf(y1[gas], Ta, Pa)
    wcvf_a['N_2'] = wcvf(1-y1[gas], Ta, Pa)
    if model == 'IAST':
      Pa_Cori = numpy.array([y1[gas], 1-y1[gas]]) * Pa
      qa['CO_2'], qa['N_2'] = pyiast.iast(
          Pa_Cori, [iso_ta['CO_2'], iso_ta['N_2']], verboseflag=False,
          #adsorbed_mole_fraction_guess=[0.99999999, 0.00000001]
      )
      logging.debug(
          "qa['CO_2']: %r, qa['N_2']: %r" % (
	    qa['CO_2'], qa['N_2']
              #, qas['CO_2']: %r, qas['N_2'], qas['CO_2'], qas['N_2']
	    )
	  )
    elif model == 'COMP':
      if not fnmatch(struc, "*dobpdc"):
        qa['CO_2'], qa['N_2'] = utils.compQ(Ta, Pa, kqv, y1[gas], dH)
      else:
        qa['CO_2'], qa['N_2'] = utils.step_load(Ta, Pa, kqv['CO_2'], y1[gas]), 0.
      logging.debug("qa['CO_2']: %r, qa['N_2']: %r" % (qa['CO_2'], qa['N_2']))
  # desorption et al.
  wcvf_d = {}
  wcvf_d['CO_2'] = wcvf(y1[1], Td, Pd)
  wcvf_d['N_2'] = wcvf(1-y1[1], Td, Pd)
  qd, qds = {}, {}
  if model == 'IAST':
    Pd_Cori = numpy.array([y1[1], 1-y1[1]]) * Pd
    qd['CO_2'], qd['N_2'] = pyiast.iast(
      Pd_Cori, [iso_td['CO_2'], iso_td['N_2']], verboseflag=False,
      adsorbed_mole_fraction_guess=[0.99999999, 0.00000001]
    )
    logging.debug(
      "qd['CO_2']: %r, qd['N_2']: %r" % (
        qd['CO_2'], qd['N_2']
        #, qds['CO_2']: %r, qds['N_2'], qds['CO_2'], qds['N_2']
      )
    )
  elif model == 'COMP':
    if not fnmatch(struc, "*dobpdc"):
      qd['CO_2'], qd['N_2'] = utils.compQ(Td, Pd, kqv, y1[1], dH)
    else:
      qd['CO_2'], qd['N_2'] = utils.step_load(Td, Pd, kqv['CO_2'], y1[1]), 0.
    logging.debug("qd['CO_2']: %r, qd['N_2']: %r" % (qd['CO_2'], qd['N_2']))
  wc, wct = {}, {} #wcs = {}
  for m in ['CO_2', 'N_2']:
    wc[m] = qa[m] - qd[m]
    wct[m] = wc[m] + wcvf_a[m] - wcvf_d[m]
  if wc['CO_2'] < 0: # or wcs['CO_2'][0] < 0 or wcs['CO_2'][1] < 0: #check validity
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
    #if model == 'IAST':
    #  logging.debug("wcs['CO_2'] = %r, wcs['N_2'] = %r" % (wcs['CO_2'], wcs['N_2']))
    return Qs - Qd, pur, m1, wct['CO_2'], Qs, -Qd

def totalW(Pd, pur):
  plog = -5.4433e5/math.log(1e6/1e3)
  Wp = 7.1723e5+math.log(Pd/1000)*plog
  mp = (2.313-2.102)/(math.sqrt(1e6)-math.sqrt(1e3))
  mpr = 1.102+mp*(math.sqrt(Pd)-math.sqrt(1e3))
  Wt = Wp*(1+mpr*(1./pur-1))
  return Wt

def main(args):
  global cp, ms, Ta, Pa, dH, vf, gas, model, kqv, totE, pCO2, method
  global iso_ta, iso_df, T_df, iso_td
  db = ccsdb.load(args.db)
  model = args.model
  props = ccsdb.getProps(db, args.struc)
  method = args.method
  limit = args.limit
  if limit == 'LIMIT_UP':
    cp, rho = float(props[3]), float(props[1])
  elif limit == 'LIMIT_LOW':
    cp, rho = float(props[2]), float(props[1])
  else:
    cp, rho = float(props[0]), float(props[1])
  V, vf = ccsdb.getHead(db, 'FixedBed')
  ms = rho * ( 1. - vf ) # V = 1m3
  logging.debug('ms = %f' %  ms)
  FlueGasHead = ccsdb.getHead(db, 'FlueGas')
  Ta, Pa = [ float(a) for a in FlueGasHead ]
  dH['CO_2'], dH['N_2'] = ccsdb.getHoAs(db, args.struc)
  if args.comp == 'coal': gas = 0; totE = 6631.2 ; pCO2 = 1.80
  elif args.comp == 'NG': gas = 2; totE = 21023.26; pCO2 = 3.23
  else:
    gas = 3
    Ta = float(288)
  for m in ccsdb.getAdsorbates(db, args.struc):
    kqv[m] = ccsdb.getAdsorbate(db, args.struc, m).get('kq300')
    df_dir = os.path.join('ccsdata', args.struc, m)
    df_file = glob(os.path.join(df_dir, '*.csv'))[0]
    T_df[m] = int(os.path.splitext(os.path.basename(df_file))[0][:-1])
    iso_df[m] = pd.read_csv(df_file, sep=' ')
    iso_ta[m], _ = deltaT.ConstructInterpolatorIsothermAtTnew(
      iso_df[m], T_df[m], Ta, loading_key="loading(mol/kg)", pressure_key="pressure(Pa)",
      hoa_key="HoA(kJ/mol)", fill_value=iso_df[m]["loading(mol/kg)"].max()
    )

  #--------
  # minimization steps:
  #Pd = [0.01:3] in 1.0 atm
  #Td = [333:473] in 10 K
  Pstep_size = 0.01
  Tstep_size = 10
  aPE = []
  cnt = 0
  Toffset = 308 if args.comp == 'air' else 333
  if limit == 'IAST_PSA':
    Td = 333
    data = []
    for m in ccsdb.getAdsorbates(db, args.struc):
      iso_td[m], _ = deltaT.ConstructInterpolatorIsothermAtTnew(
        iso_df[m], T_df[m], Td, loading_key="loading(mol/kg)", pressure_key="pressure(Pa)",
        hoa_key="HoA(kJ/mol)", fill_value=iso_df[m]["loading(mol/kg)"].max()
      )
    for p in range(30):
      Pd = (p*Pstep_size+0.01)*1013.25/0.01 #convert from atm to Pa 1013.25*10^2
      logging.debug("Td=%r, Pd=%r" % (Td, Pd))
      #if not cnt%100:
      #  sys.stdout.write('.')
      #  sys.stdout.flush()
      PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt = totalPE(db, args.struc, float(Td), Pd)
      logging.debug("PE: %r" %  PE)
      if PE < 0: continue
      data.append([Td, Pd, PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt])
      cnt += 1
    data = numpy.array(data)
    if len(data) > 0:
      optPE = data[numpy.argmin(data.T[2])]
      aPE.append(optPE)
  elif limit == 'IAST_TSA':
    Pd = 10132.50
    Td = 353
    data = []
    for m in ccsdb.getAdsorbates(db, args.struc):
      iso_td[m], _ = deltaT.ConstructInterpolatorIsothermAtTnew(
        iso_df[m], T_df[m], Td, loading_key="loading(mol/kg)", pressure_key="pressure(Pa)",
        hoa_key="HoA(kJ/mol)", fill_value=iso_df[m]["loading(mol/kg)"].max()
      )
    logging.debug("Td=%r, Pd=%r" % (Td, Pd))
    #if not cnt%100:
    #  sys.stdout.write('.')
    #  sys.stdout.flush()
    PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt = totalPE(db, args.struc, float(Td), Pd)
    logging.debug("PE: %r" %  PE)
    #if PE < 0: continue
    data.append([Td, Pd, PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt])
    cnt += 1
    data = numpy.array(data)
    if len(data) > 0:
      optPE = data[numpy.argmin(data.T[2])]
      aPE.append(optPE)
  else:
    #Td = 333
    #Pd = 101325.0
    for t in range(100): #[0]: #range(70):
      Td = t*Tstep_size + Toffset
      for m in ccsdb.getAdsorbates(db, args.struc):
        iso_td[m], _ = deltaT.ConstructInterpolatorIsothermAtTnew(
          iso_df[m], T_df[m], Td, loading_key="loading(mol/kg)", pressure_key="pressure(Pa)",
          hoa_key="HoA(kJ/mol)", fill_value=iso_df[m]["loading(mol/kg)"].max()
        )
      data = []
      for p in range(30):
        Pd = (p*Pstep_size+0.01)*1013.25/0.01 #convert from atm to Pa 1013.25*10^2
        logging.debug("Td=%r, Pd=%r" % (Td, Pd))
        if not cnt%100:
          sys.stdout.write('.')
          sys.stdout.flush()
        PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt = totalPE(db, args.struc, float(Td), Pd)
        logging.debug("PE: %r" %  PE)
        if PE < 0: continue
        data.append([Td, Pd, PE, Cap, Comp, MProd, wct, Qs, Qd, pur, Qt])
        cnt += 1
      data = numpy.array(data)
      if len(data) > 0:
        optPE = data[numpy.argmin(data.T[2])]
        aPE.append(optPE)
  aPE = numpy.array(aPE)
  logging.debug("aPE: %r" %  aPE)
  optPET = aPE[numpy.argmin(aPE.T[2])].tolist()
  finPE = utils.trimNumber(optPET[2]/1000., 6)
  finCap = utils.trimNumber(optPET[3]/1000., 6)
  finComp = utils.trimNumber(optPET[4]/1000., 6)
  finQs = utils.trimNumber(optPET[7]/1000., 6)
  finQd = utils.trimNumber(optPET[8]/1000., 6)
  finMProd = utils.trimNumber(optPET[5], 5)
  finWC = utils.trimNumber(optPET[6], 5)
  finPur = utils.trimNumber(optPET[9], 5)
  finQt = utils.trimNumber(optPET[10]/1000., 6)
  finP = utils.trimNumber(optPET[1], 5)
  finT = optPET[0]
  if (args.comp == 'coal') or (args.comp == 'NG'):
    finEL = utils.trimNumber(finPE * pCO2 / totE, 5)
  else:
    finEL = -1
  #--------

  finPE_id = [finPE, finP, finT, finEL, finCap, finComp, finMProd, finWC, finPur]
  logging.debug("%s %r" %  (args.struc, finPE_id))
  print (args.struc, finPE_id)
  ccsdb.writePEs(db, str(args.struc), str(args.comp), str(args.limit), finPE_id)
  ccsdb.dump(db, args.db)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-l", "--log", help="create log file", action = "store_true")
  parser.add_argument("db", help="db.yml")
  parser.add_argument("struc", help="CuBTC")
  parser.add_argument("comp", help="coal, NG or air")
  parser.add_argument("model", help="IAST or COMP")
  parser.add_argument("limit", help="IAST, COMP, IAST_PSA, IAST_TSA, LIMIT_UP, LIMIT_LOW")
  parser.add_argument("method", nargs='?', default="carnot", help="linear")
  args = parser.parse_args()
  if args.log:
    now = datetime.now().isoformat()
    logfile = 'calPE_%s_%s.log' % (args.struc, now)
    logging.basicConfig(
      filename=logfile, format='%(message)s', level=logging.DEBUG
    )
  main(args)
