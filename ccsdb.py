#!/usr/bin/env python

import fnmatch
import os
import yaml
import numpy

# ordered dict for yaml
# ---------------------
# http://stackoverflow.com/questions/8651095/controlling-yaml-serialization-order-in-python
# from collections import OrderedDict
# http://stackoverflow.com/questions/14358162/funnelweb-error-cannot-import-ordereddict
try:
  from collections import OrderedDict
except ImportError:
  # python 2.6 or earlier, use backport
  from ordereddict import OrderedDict

def dump_as_map(adict):
  yaml.add_representer(adict, _represent_dictorder)

def _represent_dictorder(self, data):
  if isinstance(data, Structure) or isinstance(data, Adsorbate) or isinstance(data, Database) or isinstance(data, ParaEnergy) or isinstance(data, IDpe) or isinstance(data, EvalCrit) or isinstance(data, DbHeader):
    return self.represent_mapping('tag:yaml.org,2002:map', data.__getstate__().items())
  else:
    return self.represent_mapping('tag:yaml.org,2002:map', data.items())

# yaml classes
# ------------
class Adsorbate(yaml.YAMLObject):
  yaml_tag = u'!Adsorbate'
  def __init__(self,T,kH1,qs1,v1,kH2,qs2,v2,dH,kq):
    self.T = T
    self.kH1 = kH1
    self.qs1 = qs1
    self.v1 = v1
    self.kH2 = kH2
    self.qs2 = qs2
    self.v2 = v2
    self.dH = dH
    self.kq300 = kq
  def __getstate__(self):
    d = OrderedDict()
    d['T'] = self.T
    d['kH1'] = self.kH1
    d['qs1'] = self.qs1
    d['v1'] = self.v1
    d['kH2'] = self.kH2
    d['qs2'] = self.qs2
    d['v2'] = self.v2
    d['dH'] = self.dH
    d['kq300'] = self.kq300
    return d
  def __repr__(self):
    return "%s(T=%r,kH1=%r,qs1=%r,v1=%r,kH2=%r,qs2=%r,v2=%r,dH=%r,kq300=%r)" % (
      self.__class__.__name__, self.T,
      self.kH1, self.qs1, self.v1, self.kH2, self.qs2, self.v2, self.dH, self.kq300
    )

class DbHeader(yaml.YAMLObject):
  yaml_tag = u'!DbHeader'
  def __init__(self,datadir,bed,cond,comp,line):
    self.datadir = datadir
    self.bed = bed
    self.cond = cond
    self.comp = comp
    self.line = line
  def __getstate__(self):
    d = OrderedDict()
    d['DataDir'] = self.datadir
    d['FixedBed'] = self.bed
    d['FlueGas'] = self.cond
    d['IastY1'] = self.comp
    d['EnvelopeLine'] = self.line
    return d
  def __repr__(self):
    return "%s(datadir=%r,bed=%r,cond=%r,comp=%r,line=%r)" % (
	self.__class__.__name__, self.datadir, self.bed,
	self.cond, self.comp, self.line
	)

class ParaEnergy(yaml.YAMLObject):
  yaml_tag = u'!ParaEnergy'
  def __init__(self,coal,NG,air):
    self.coal = coal
    self.NG = NG
    self.air = air
  def __getstate__(self):
    d = OrderedDict()
    d['coal'] = self.coal
    d['NG'] = self.NG
    d['air'] = self.air
    return d
  def __repr__(self):
    return "%s(coal=%r,NG=%r,air=%r)" % (
	self.__class__.__name__, self.coal, self.NG, self.air)

class IDpe(yaml.YAMLObject):
  yaml_tag = u'!ParaEnergyIDs'
  def __init__(self,COMP,IAST,IAST_PSA,IAST_TSA,LIMIT_UP, LIMIT_LOW):
    self.COMP = COMP
    self.IAST = IAST
    self.IAST_PSA = IAST_PSA
    self.IAST_TSA = IAST_TSA
    self.LIMIT_UP = LIMIT_UP
    self.LIMIT_LOW = LIMIT_LOW
  def __getstate__(self):
    d = OrderedDict()
    d['COMP'] = self.COMP
    d['IAST'] = self.IAST
    d['IAST_PSA'] = self.IAST_PSA
    d['IAST_TSA'] = self.IAST_TSA
    d['LIMIT_UP'] = self.LIMIT_UP
    d['LIMIT_LOW'] = self.LIMIT_LOW
    return d
  def __repr__(self):
    return "%s(COMP=%r,IAST=%r,IAST_PSA=%r,IAST_TSA=%r,LIMIT_UP=%r,LIMIT_LOW=%r)" % (
	self.__class__.__name__, self.COMP, self.IAST, self.IAST_PSA, self.IAST_TSA, self.LIMIT_UP, self.LIMIT_LOW)

class EvalCrit(yaml.YAMLObject):
  yaml_tag = u'!EvalCrit'
  def __init__(self,up,sel,wc,reg,sorb,pe):
    self.Uptake = up
    self.Selectivity = sel
    self.WorkCap = wc
    self.Regen = reg
    self.SorbSel = sorb
    self.PE = pe
  def __getstate__(self):
    d = OrderedDict()
    d['Uptake'] = self.Uptake
    d['Selectivity'] = self.Selectivity
    d['WorkCap'] = self.WorkCap
    d['Regenerate'] = self.Regen
    d['SorbSel'] = self.SorbSel
    d['PE'] = self.PE
    return d
  def __repr__(self):
    return "%s(Uptake=%r,Selectivity=%r,WorkCap=%r,Regenerate=%r,SorbSel=%r,PE=%r)" % (
	self.__class__.__name__, self.Uptake, self.Selectivity, self.WorkCap, self.Regen, self.SorbSel, self.PE)

class Structure(yaml.YAMLObject):
  yaml_tag = u'!Structure'
  def __init__(self,pe,A,p,ec):
    self.PE = pe
    self.Adsorbates = A
    self.Props = p
    self.EvalCrit = ec
  def __getstate__(self):
    d = OrderedDict()
    d['Adsorbates'] = self.Adsorbates
    d['PE'] = self.PE
    d['Props'] = self.Props
    d['EvalCrit'] = self.EvalCrit
    return d
  def __repr__(self):
    return "%s(PE=%r,A=%r,Props=%r,EvalCrit=%r)" % (
      self.__class__.__name__, self.PE, self.Adsorbates, self.Props, self.EvalCrit)

class Database(yaml.YAMLObject):
  yaml_tag = u'!Database'
  def __init__(self,name,head,S):
    self.Header = head
    self.Structures = S
  def __getstate__(self):
    d = OrderedDict()
    d['Header'] = self.Header
    d['Structures'] = self.Structures
    return d
  def __repr__(self):
    return "%s(head=%r,S=%r)" % (
      self.__class__.__name__, self.Header,
      self.Structures)

# detect differences in dictionary
# --------------------------------
class DictDiffer(object):
  """
  Calculate the difference between two dictionaries as:
  (1) items added
  (2) items removed
  (3) keys same in both but changed values
  (4) keys same in both and unchanged values
  """
  def __init__(self, current_dict, past_dict):
    self.current_dict, self.past_dict = current_dict, past_dict
    self.current_keys, self.past_keys = [set(d.keys()) for d in (current_dict, past_dict)]
    self.intersect = self.current_keys.intersection(self.past_keys)
  def added(self):
    return self.current_keys - self.intersect
  def removed(self):
    return self.past_keys - self.intersect
  def changed(self):
    return set(o for o in self.intersect if self.past_dict[o] != self.current_dict[o])
  def unchanged(self):
    return set(o for o in self.intersect if self.past_dict[o] == self.current_dict[o])

# database
# -------
# load database
def load(f):
  stream = open(f, 'r')
  return yaml.load(stream)

# dump the database to file
def dump(db, f):
  dumpMaps()
  stream = open(f, 'w')
  return yaml.dump(db, stream, explicit_start=True, indent=2)

# access database content
# -----------------------
# get list of header (dict)
def getHeader(db):
  return db['Header']

# get header (dictionary)
def getHead(db, head):
  return getHeader(db).get(head)

# get list of structures (dict)
def getStructures(db):
  return db['Structures']

# get structure (dictionary)
def getStructure(db, struc):
  return getStructures(db).get(struc)

# get list of PEs for gas compositions (dict)
def getPEs(db, struc):
  return getStructure(db, struc).get('PE')

# get PEs for gas composition (dictionary)
def getPE(db, struc, gas):
  return getPEs(db, struc).get(gas)

# get PE for specific ID (array: PE, P, T)
def getPE_id(db, struc, gas, pe_id):
  return getPE(db, struc, gas).get(pe_id)

# get list of EvalCrits for structure (dict)
def getEvalCrits(db, struc):
  return getStructure(db, struc).get('EvalCrit')

# get EvalCrits (dictionary)
def getEvalCrit(db, struc, crit):
  return getEvalCrits(db, struc).get(crit)

# get list of adsorbates (dict)
def getAdsorbates(db, struc):
  return getStructure(db, struc).get('Adsorbates')

# get adsorbate (dictionary)
def getAdsorbate(db, struc, mol):
  return getAdsorbates(db, struc).get(mol)

# get temperatures (array)
def getT(db, struc, mol):
  return getAdsorbate(db, struc, mol).get('T')

# set temperatures new for update
def setT(db, struc, mol, T):
  a = getAdsorbate(db, struc, mol)
  a['T'] = T

# get file input path (*.csv)
def getInputPath(db, struc, mol, t):
  fn = str(t) + 'K.csv'
  return os.path.join(getHead(db, 'DataDir'), struc, mol, fn)

# write results of isotherm fits
def writeIsoFit(db, struc, mol, kH1, qs1, v1, kH2, qs2, v2, dH, kq300):
  ads = getAdsorbate(db, struc, mol)
  ads['kH1'] = kH1
  ads['qs1'] = qs1
  ads['v1'] = v1
  ads['kH2'] = kH2
  ads['qs2'] = qs2
  ads['v2'] = v2
  ads['dH'] = dH
  ads['kq300'] = kq300
  return ads

# write results of calPE
def writePEs(db, struc, gas, pe_id, finPE_id):
  pe = getPE(db, struc, gas)
  pe[pe_id] = finPE_id
  return pe

# write results of EnvelopeLine
def writeEnvelope(db, gas, pars):
  env = getHead(db, 'EnvelopeLine')
  env[gas] = pars
  return env

# write results of evalcrit
def writeEvalCrits(db, struc, q, sel, wc, reg, sorb, pe):
  evalcrit = getEvalCrits(db, struc)
  evalcrit['Uptake'] = [float(e) for e in q]
  evalcrit['Selectivity'] = [float(s) for s in sel]
  evalcrit['WorkCap'] = [float(w) for w in wc]
  evalcrit['Regenerate'] = float(reg)
  evalcrit['SorbSel'] = float(sorb)
  evalcrit['PE'] = float(pe)
  return evalcrit

# get HoA's
def getHoAs(db, struc):
  return getAdsorbate(db, struc, 'CO_2').get('dH'), getAdsorbate(db, struc,
                                                                 'N_2').get('dH')#, getAdsorbate(db, struc, 'H2O').get('dH')#, getAdsorbate(db, struc, 'C3H8').get('dH'), getAdsorbate(db, struc, 'H2S').get('dH')#, getAdsorbate(db, struc, 'H2O').get('dH')

# load data file
def loadData(db, struc, mol, t):
  return numpy.loadtxt(getInputPath(db, struc, mol, t))

# get Structure Props (cp, rho)
# -----------------------------
def getProps(db, struc):
  return getStructure(db, struc).get('Props')

def dumpMaps():
  dump_as_map(Database);
  dump_as_map(Structure);
  dump_as_map(Adsorbate);
  dump_as_map(ParaEnergy);
  dump_as_map(IDpe);
  dump_as_map(EvalCrit);
  dump_as_map(DbHeader);

# update database (scans the datdir)
# --------------------------------------------------
def update(f, datdir):
  datdir = os.path.realpath(datdir)
  if not os.path.exists(f):
    fg = [313., 1.01325e5]
    fb = [40.0, 0.35]
    y1 = [0.14, 0.99, 0.04, 0.0004]
    db = load('../ccsscripts/NEWoptPEdb.yml')
    env = getHead(db, 'EnvelopeLine')
    head = DbHeader(datdir, fb, fg, y1, env)
    S = {}
  else:
    db = load(f) # else load f
    #env = ParaEnergy([-1]*4,[-1]*4,[-1]*4)
    #head = getHeader(db)
    #head['EnvelopeLine'] = env
    S = getStructures(db)
  for struc in os.listdir(datdir): #getStructures(db):
    if fnmatch.fnmatch(struc, '.git'): continue
    if fnmatch.fnmatch(struc, '.txt-files'): continue
    if struc in S: continue # attention: comment out if new pars
    curStruc = os.path.join(datdir,struc)
    if os.path.isdir( curStruc ):
      if not struc in S: # init Props, PEs and EvalCrits with defaults
        pr = [0.985e3, 1.0, 0.761e3, 1.210e3]
        pe_coal = IDpe([-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8)
        pe_ng = IDpe([-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8)
        pe_air = IDpe([-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8)
        PE = ParaEnergy(pe_coal,pe_ng,pe_air)
        up = [-1, -1]
        sel = [-1, -1]
        wc = [-1, -1]
        reg = [-1]
        sorb = [-1]
        pe = [-1]
        EC = EvalCrit(up,sel,wc,reg,sorb,pe)
      else: # copy Props, PEs and EvalCrits
        pr = getProps(db, struc)
        pe_air = getPEs(db, struc).get('air')
        #pe_airNew = IDpe(pe_air['COMP'], pe_air['IAST'], pe_air['IAST_PSA'], pe_air['IAST_TSA'], [-1]*5, [-1]*5)
        pe_coal = getPEs(db, struc).get('coal')
        #pe_coalNew = IDpe(pe_coal['COMP'], pe_coal['IAST'], pe_coal['IAST_PSA'], pe_coal['IAST_TSA'], [-1]*5, [-1]*5)
        pe_ng = getPEs(db, struc).get('NG')
        #pe_water = IDpe([-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8,[-1]*8)
        #pe_ngNew = IDpe(pe_ng['COMP'], pe_ng['IAST'], pe_ng['IAST_PSA'], pe_ng['IAST_TSA'], [-1]*5, [-1]*5)
        PE = ParaEnergy(pe_coal,pe_ng,pe_air)
        EC = getEvalCrits(db,struc)
      A = {} # adsorbates
      for mol in ['CO_2', 'N_2']:#, 'H2O']: #os.listdir( curStruc ):
        curMol = os.path.join(curStruc, mol)
        T = [] #temperatures (always keep all temperatures!)
        for temp in os.listdir(curMol):
          if fnmatch.fnmatch(temp,'*.csv'):
            T.append( int(temp[:-5]) )
        if not struc in S:
          A[mol] = Adsorbate(T,[],[],[],[],[],[],[],[])
        else:
          oldA = getAdsorbate(db, struc, mol)
          if oldA is None: A[mol] = Adsorbate(T,[],[],[],[],[],[],[],[])#continue
          else:
            T = oldA.get('T')
            kH1 = oldA.get('kH1')
            qs1 = oldA.get('qs1')
            v1 = oldA.get('v1')
            kH2 = oldA.get('kH2')
            qs2 = oldA.get('qs2')
            v2 = oldA.get('v2')
            dH = oldA.get('dH')
            kq = oldA.get('kq300')
            A[mol] = Adsorbate(T,kH1,qs1,v1,kH2,qs2,v2,dH,kq)
      S[struc] = Structure(PE,A,pr,EC)
  if not os.path.exists(f):
    db = Database(f,head,S)
  dump(db,f)
  #print yaml.dump(Database(f,datdir,fg,fb,y1,S), explicit_start=True, indent=2)

# main function for execution as script
# -------------------------------------
if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("dbfile", help="path to db file")
  parser.add_argument("datdir", help="path to data directory")
  args = parser.parse_args()
  update(args.dbfile, args.datdir)
