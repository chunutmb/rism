#!/usr/bin/env python


import os, sys, re, getopt


rmin = 4
rmax = 16
rdel = 0.1
fncfg = "NaCl.cfg"
progdir = "../../prog"
prog0 = "rism0"
prog = None


def help():
  ''' print the help message and quit '''
  prog = sys.argv[0]
  print "%s [Options] template.cfg\n" % prog
  print " Options:"
  print "   --rmin:     followed by the minimal radius"
  print "   --rmax:     followed by the maximal radius"
  print "   --rdel:     followed by the radius increment"
  print "   --progdir:  followed by the directory that contains", prog0
  exit(1)



def doargs():
  global rmin, rmax, rdel, fncfg
  global progdir, prog

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "h",
        [ "rmin=", "rmax=", "rdel=", "progdir=",
          "help" ] )
  except getopt.GetoptError, err:
    help()

  for o, a in opts:
    if o in ("--rmin",):
      rmin = float(a)
    elif o in ("--rmax",):
      rmax = float(a)
    elif o in ("--rdel",):
      rdel = float(a)
    elif o in ("--progdir",):
      progdir = a
    elif o in ("-h", "--help"):
      help()

  if len(args) > 0:
    fncfg = args[0]

  prog = os.path.join(progdir, prog0)



def loadcfg(fn):
  ''' load the configuration template '''
  s = open(fn).readlines()
  n = len(s)
  # search for the last line that specifies the distance
  for i in range(n-1, -1, -1):
    ln = s[i].strip()
    if ln.startswith("dis("): break
  if i < 0:
    print "no distance matrix line found in", fn
    raise Exception
  return s, i



def runscript(cfg0, cfgln, r, name = None):
  ''' change the distance of the special pair to r
      then run the script
      return the chemical potential '''

  print "running %s with r = %s..." % (prog0, r)

  # construct an input file
  cfg = cfg0[:]
  if r == None:
    cfg[cfgln] = "\n"
    opts = ""
  else:
    arr = cfg[cfgln].split("=") # parse it into `key = value'
    arr[1] = (" %s" % r) # change the value part
    cfg[cfgln] = "=".join(arr)
    opts = "-!" # this option disables the solute-solute interaction
  fncfg = "test_r%s.cfg" % r
  fnout = "out_r%s.dat" % r
  open(fncfg, "w").writelines(cfg)

  cmd = "%s %s %s > %s" % (prog, opts, fncfg, fnout)

  # run the command
  ret = os.system(cmd)
  if ret != 0:
    print "command [%s] failed with %s" % (cmd, ret)
    exit(1)

  # analyze the output
  bmu = 0
  if r != None:
    s = open(fnout).readlines()
    n = len(s)
    for i in range(n - 1, -1, -1):
      ln = s[i]
      m = re.search("mol.*chem.*pot.*X beta(.*)", ln)
      if m:
        bmu = float( m.group(1).strip() )
        break

  if r == None:
    os.rename("out.dat", name + "_out.dat")
  os.remove(fnout)
  os.remove(fncfg)
  os.remove("crdnum.dat")
  return bmu



def doit():
  # build the program
  os.system("make -C " + progdir)

  # load the configuration file
  cfgtemplate, cfgln = loadcfg(fncfg)
  name = fncfg.split(".")[0]

  i = int( rmin / rdel );
  r = i * rdel
  ls = []
  while r < rmax + .5 * rdel:
    # change the distance to r and run the script
    bmu = runscript(cfgtemplate, cfgln, r)
    ls += [(r, bmu), ]
    i += 1
    r = i * rdel

  # run the script with infinite separation
  bmu0 = runscript(cfgtemplate, cfgln, 1e30)
  sout = "#    r    beta*mu\n"
  for i in range(len(ls)):
    ls[i] = (ls[i][0], ls[i][1] - bmu0)
    sout += "%9s %s\n" % (ls[i][0], ls[i][1])

  # run the script without the distance constraint
  runscript(cfgtemplate, cfgln, None, name)

  fnbmu = fncfg.split(".")[0] + "_bmu.dat"
  open(fnbmu, "w").write(sout)



if __name__ == "__main__":
  doargs()
  doit()
