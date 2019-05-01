#!/usr/bin/python3.4
 
import QE as qe
import sys
import numpy as np
 
class Files:
  def __init__(self):
    self.prefix = ''
    self.pdos = """&projwfc
    Emin=-100, Emax=100.0, DeltaE=0.1
    ngauss=1, degauss=0.02
 /
"""
    self.dos = """&dos
  Emin=-100.0, Emax=100.0, DeltaE=0.1
/
"""
    self.fermi = """&fermi
  Emin=-100.0, Emax=100.0, DeltaE=0.1
/
"""
    self.bands = """&bands
  Emin=-100.0, Emax=100.0, DeltaE=0.1
/
"""


    self.scf = """"""
    self.nscf = """"""
    self.round = 5
  def Process(self,filestring,QE):
    try:
      f = open(filestring,'r')
    except:
      print("Cannot open %s" % filestring)
      sys.exit(1)
    self.scf = ''
    self.nscf = ''
    for i in f:
      if 'calculation' in i:
        self.scf += i.replace('vc-relax','scf').replace('relax','scf')
        self.nscf += i.replace('vc-relax','nscf').replace('relax','nscf')
      elif '&system' in i:
        tmpstr = '''&system\n  ibrav = 14\n'''
        for j in ['a','b','c']:
          tmpstr += "  " + j + " = " + str(round(QE.norms[j],self.round)) + "\n"
        tmpstr += "  cosBC = " + str(round(np.cos(QE.angles['alpha'] * np.pi/180.),self.round)) + "\n"
        tmpstr += "  cosAC = " + str(round(np.cos(QE.angles['beta'] * np.pi/180.),self.round)) + "\n"
        tmpstr += "  cosAB = " + str(round(np.cos(QE.angles['gamma'] * np.pi/180.),self.round)) + "\n"
        self.scf += tmpstr
        self.nscf += tmpstr
        for j in range(0,7):
          next(f)
      elif 'ATOMIC_POSITIONS' in i.upper():
        tmpstr = 'ATOMIC_POSITIONS {angstrom}\n'
        for atom, loc  in QE.atoms.items():
          tmpat = ''.join(i for i in atom if not i.isdigit())
          tmpstr += tmpat + "\t" + str(loc[0]) + "\t" + str(loc[1]) + "\t" + str(loc[2]) + "\n"
          next(f)
        self.scf += tmpstr
        self.nscf += tmpstr
      elif 'K_POINTS' in i.upper():
        self.scf += i 
        self.nscf += i 
        line = next(f)
        self.scf += line
        ks = line.split()
        for j in range(0,6):
          if j < 3:
            self.nscf += str(24) + " "
          else:
            self.nscf += ks[j] + " "
        self.nscf += "\n"
      elif "prefix" in i.lower():
        self.prefix=i.replace('prefix','').replace('=','').replace("'","").replace('"','').strip()
      else:
        self.scf += i
        self.nscf += i
      self.pdos = self.pdos.replace('PREFIX',self.prefix)
      self.dos = self.dos.replace('PREFIX',self.prefix)
         
def main(command):
  QEStruct= qe.Struct()
  QEStruct.File_Process(command[1])
  outputs = Files()
  outputs.Process(command[0],QEStruct)
  f = open('scf.in','w')
  f.write(outputs.scf)
  f.close()
  f = open('nscf.in','w')
  f.write(outputs.nscf)
  f.close()
  f = open('dos.in','w')
  f.write(outputs.dos)
  f.close()
  f = open('pdos.in','w')
  f.write(outputs.pdos)
  f = open('fermi.in', 'w')
  f.write(outputs.fermi)
  f.close()
 
if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Incorrect number of arguments, run as ./Scf.py QEINPUT QEOUTPUT")
    sys.exit(6)
  command = [sys.argv[1],sys.argv[2]]
  main(command)
