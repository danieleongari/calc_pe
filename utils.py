#!/usr/bin/env python

import os
import math, copy
import numpy as np
import pandas as pd
import pyiast

R = 8.3144621 # J/K/mol

colarr = [ "0", "1", "2", "3", "4", "5", "6", "9",
           "rgb \"#ff8c00\"", "rgb \"#228b22\"" , "rgb \"#b22222\"",
           "rgb \"#9370db\"", "rgb \"#bdb76b\"", "rgb \"#00bfff\"",
           "rgb \"#fa8072\"", "rgb \"#ee82ee\"", "rgb \"#7fffd4\"",
           "rgb \"#0000cd\"", "rgb \"#ffdab9\"", "rgb \"#eee9e9\"",
           "rgb \"#eecbad\"", "rgb \"#a0522d\"", "rgb \"#2e8b57\"",
           "rgb \"#3cb371\"", "rgb \"#20b2aa\"", "rgb \"#98fb98\"",
           "rgb \"#db7093\"", "rgb \"#b03060\"", "rgb \"#c71585\"",
           "rgb \"#bc8f8f\"", "rgb \"#cd5c5c\"" ]

globkey = [
  'spacing 1.4', 'samplen 1.2', 'reverse Left',
  'box lw 2', 'height 0.5'
]

def mypow(x, a, b): return a*np.power(x,b)

def mypow2(x, p, s): # x = press, p = par list for s = site index
  return mypow(x, p[s*3] / p[s*3+1], p[s*3+2])

# competitive langmuir isotherm
# kqv['CO_2'][s*3+p]: mol & site index (s) &
# parameter index (p) [ (k,0) (q,1) (v,2) ]
def compQ(T, P, p, y, h):
  pdc = copy.deepcopy(p)
  for m in ['CO_2', 'N_2']:
    for s in [0,1]:
      pdc[m][s*3] = calckH(T, p[m][s*3], h[m][s])
  isDual = (p['CO_2'][3] > 0 and p['CO_2'][4] > 0)
  Q = {}
  for m in ['CO_2', 'N_2']:
    Qdenom = 1
    if not isDual:
      Qdenom += mypow2(P*y, pdc['CO_2'], 0)
      Qdenom += P*(1-y) * (pdc['N_2'][0]/pdc['CO_2'][1])
      #Qdenom += mypow2(P*(1-y), pdc['N_2'], 0)
    else:
      if m == 'CO_2':
        Qdenom += mypow2(P*y, pdc['CO_2'], 0)
      else:
        Qdenom += mypow2(P*y, pdc['CO_2'], 1)
        Qdenom += P*(1-y) * (pdc['N_2'][0]/pdc['CO_2'][4])
      #Qdenom += mypow2(P*(1-y), pdc['N_2'], 0)
    Q[m] = mypow(
      P*y if m == 'CO_2' else P*(1-y),
      pdc[m][0], pdc[m][2]
    )
    Q[m] /= Qdenom
    if isDual and m == 'CO_2':
      Qnumer = mypow2(P*y, pdc['CO_2'], 1)
      Qdenom = 1. + mypow2(P*y, pdc['CO_2'], 1)
      Qdenom += P*(1-y) * (pdc['N_2'][0]/pdc['CO_2'][4])
      #Qdenom += mypow2(P*(1-y), pdc['N_2'], 0)
      Q[m] += Qnumer/Qdenom
      #print "Qdenom %s %f %r" % (m, Qdenom, Q[m])
  return Q['CO_2'], Q['N_2']

# isotherm loading for step function
# (x<x0) ? a*(x/b)**c : r*(x/b)**s
# ln x0 = -m/T + t (m from dH)
def step_load(T, P, p, y):
  pdc = copy.deepcopy(p) #pars: a,c,x0_at_fit,r,s,t,h
  m = pdc[6]*1000./R
  x0 = math.exp(linear(1./T, m, pdc[5]))
  b = x0/pdc[2]
  return mypow(P*y/b, pdc[0], pdc[1]) if P*y < x0 else mypow(P*y/b, pdc[3], pdc[4])

def ext_step_load(T, P, p, y, h):
  pdc = copy.deepcopy(p) #pars: a,c,x0_at_fit,r,s,x0_at_300,Delta_b,temp
  x0 = calckH(T, pdc[5], pdc[6])
  b = pdc[2]
  return mypow(P*y/b, pdc[0], pdc[1]) if P*y < x0 else mypow(P*y/b, pdc[3], pdc[4])

# single langmuir isotherm
def single(x, k, q, v):
  return mypow(x, k, v) / ( 1 + mypow(x, k/q, v) )

# dual langmuir isotherm
def langmuir(x, k, q, v, l, r, s):
  result = single(x, k, q, v)
  if l > 0 and r > 0: result += single(x, l, r, s)
  return result

# linear interpolation in henry regime
def linear(x, m, t):
  return m*x+t

# calculate henry coefficient @T based on kH,q@300
def calckH(T, k, h):
  if not k > 0: return 0
  m = -h*1000./R
  t = math.log(k) - m/300.
  return math.exp(linear(1./T, m, t))

def linpars(x, y):
  m0 = (y[1]-y[0])/(x[1]-x[0])
  t0 = y[0]-m0*x[0]
  return m0, t0

def getEPSpath(struc, mol, t = -1):
  dir = os.path.join('figs', struc, mol)
  if not os.path.exists(dir):
    os.makedirs(dir)
  filename = struc + '_' + mol + '_'
  if t > 0:
    filename += str(t) + 'K'
  else:
    filename += 'DeltaH'
  filename += '.eps'
  return os.path.join(dir, filename)

# get step pressure for step-fcn
def getPjump(f, t):
  #temp = str(t) + 'K.csv'
  #path = str(struc) + '/' + str(mol) + '/'
  #totPath = str(struc) + '/' + str(mol) + '/' + temp
  #f = np.loadtxt('../ccsdata/' + totPath)
  global jumpPs
  jumpPs = []
  n = 0
  for n in range(len(f)-1):
    Dx = f[n+1][0] - f[n][0]
    Dy = f[n+1][1] - f[n][1]
    Pstar = (f[n][0] + f[n+1][0]) / 2.
    lnPstar = math.log(Pstar)
    #print Dx, Dy
    #return
    m = Dy / Dx
    m = math.fabs(m)
    jumpPs.append([n, m, lnPstar, Pstar])
  jumpPs = np.array(jumpPs)
  if len(jumpPs) > 0:
    maxGrad = jumpPs[np.argmax(jumpPs.T[1])]
  return t, maxGrad[2], maxGrad[3]

def ConstructInterpolatorIsothermAtTnew(df, T0, Tnew,
                                        loading_key=None,
                                        pressure_key=None,
                                        hoa_key=None,
                                        fill_value=None):
    """
    Returns InterpolatorIsotherm from pure-component isotherm data measured at
    temperature T0 extrapolated to temperature Tnew using the heat of adsorption.

    :param df: DataFrame Pandas DataFrame with pure-component isotherm tabular
               data measured at temperature T0
    :param T0: Float temperature at which isotherm in df was measured (T0)
    :param Tf: Float temperature at which we desire to extrapolate the isotherm
               in df
    :param loading_key: String key for loading column in df
    :param pressure_key: String key for pressure column in df
    :param hoa_key: String key for heat of adsorption column in df
    :param fill_value: Float value of loading to assume when an attempt is made
                       to interpolate at a pressure greater than the largest
                       pressure observed in the data

    HOA needs to be in units kJ/mol.

    Author: Cory M. Simon
    """
    if loading_key == None or pressure_key == None or hoa_key == None:
        raise Exception("Pass loading_key, hoa_key, and pressure_key," +
                        " names of loading, heat of adsorption," +
                        " and pressure cols in DataFrame.")
    if pressure_key == 'new_P':
        raise Exception("Change pressure column to something new")

    # for every point, shift pressures according to Classius-Clapyeron eqn
    R = 8.314 / 1000.0 # kJ/mol-K
    n = df.shape[0]
    df_new = pd.DataFrame()
    df_new[pressure_key] = np.zeros((n,))
    df_new[loading_key] = df[loading_key].values

    for i in range(n):
        df_new[pressure_key].iloc[i] = df[pressure_key].iloc[i] * \
                                       np.exp(-df[hoa_key].iloc[i] / R * \
                                       (1.0 / Tnew - 1.0 / T0))
    return pyiast.InterpolatorIsotherm(df_new,
                                       loading_key=loading_key,
                                       pressure_key=pressure_key,
                                       fill_value=fill_value), df_new
