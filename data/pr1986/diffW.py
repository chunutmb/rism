#!/usr/bin/env python



''' differentiate two PMF's '''


import os, sys


atom1 = 3
atom2 = 4

stem = "pr1986NaCl"
dT = 5 # temperature difference for numerical differentiation

kB = 0.00198720414667



def geteps(t):
  t -= 273.15
  eps = 87.740-0.40008*t+9.398e-4*t*t-1.410e-6*t*t*t
  deps = -0.4008+2*9.398e-4*t-3*1.410e-6*t*t
  return eps, deps



def getW(fn):
  ''' compute the potential of mean force '''
  s = open(fn).readlines()
  n = len(s)
  ri = []
  W = []
  temp = float( s[0].strip()[1:].split()[0] )
  eps, deps = geteps(temp)
  kT = temp*kB
  for i in range(n):
    ln = s[i].strip()
    if not ln or ln.startswith("#"): continue
    arr = ln.split()
    if len(arr) < 10: continue
    a1 = int(arr[6])
    a2 = int(arr[7])
    if a1 != atom1 or a2 != atom2: continue
    ri += [float(arr[0]),]
    dW = float(arr[9])
    bu = float(arr[5])
    buc = float(arr[8])
    W += [(dW + bu - buc + buc/eps)*kT,]
  return ri, W, temp



def diff(fn1, fn2):
  ''' fn2 - fn1 '''
  r1, W1, T1 = getW(fn1)
  r2, W2, T2 = getW(fn2)
  n = len(W1)
  Tav = (T1 + T2)/2
  dT = T2 - T1
  nTS = [0]*n
  E = [0]*n
  Wav = [0]*n
  s = ""
  for i in range(n):
    nTS[i] = Tav*(W2[i] - W1[i])/dT
    Wav[i] = (W1[i] + W2[i])/2
    E[i] = Wav[i] - nTS[i]
    s += "%g %g %g %g %g %g\n" % (r1[i], Wav[i], nTS[i], E[i], W1[i], W2[i])
  return s



def getoutput(stem, tag, dT):
  fncfg0 = stem + ".cfg"
  fncfg = stem + tag + ".cfg"
  s = open(fncfg0).readlines()
  for i in range(len(s)):
    ln = s[i]
    arr = ln.strip().split("=")
    if arr[0].strip() == "T":
      temp = float(arr[1].strip()) + dT
      s[i] = arr[0] + "= " + str(temp) + "\n"
      break
  open(fncfg, "w").writelines(s)
  print "generating", fncfg
  fnout = stem + tag + "_out.dat"
  cmd = "../../prog/rism0 %s -o %s" % (fncfg, fnout)
  print "running", cmd
  os.system(cmd)
  return fnout



def autodiff(stem, dT):
  fnout1 = getoutput(stem, "T1", -dT)
  fnout2 = getoutput(stem, "T2", +dT)
  s = diff(fnout1, fnout2)
  fndiff = stem + "_diff.dat"
  print "writing", fndiff
  open(fndiff, "w").write(s)
  return s



if __name__ == "__main__":
  if len(sys.argv) >= 2:
    stem = sys.argv[1]
  if len(sys.argv) >= 3:
    dT = float(sys.argv[2])
  autodiff(stem, dT)


