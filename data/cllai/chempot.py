#!/usr/bin/env python



import os, sys, re, getopt, glob



rmin = 4
rmax = 16
rdel = 0.5
fncfg = "my_ions.cfg"
progdir = None
prog0 = "rism0"
prog = None
epsv = 78


def help():
  ''' print the help message and quit '''
  prog = sys.argv[0]
  print "%s [Options] template.cfg\n" % prog
  print " Options:"
  print "   --rmin=:          followed by the minimal radius"
  print "   --rmax=:          followed by the maximal radius"
  print "   --dr=, --rdel=:   followed by the radius increment"
  print "   --eps=:           followed by the solvent dielectric contant, default", epsv
  print "   --progdir=:       followed by the directory that contains", prog0
  exit(1)



def doargs():
  global rmin, rmax, rdel, fncfg, progdir, prog

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "h",
        [ "rmin=", "rmax=", "rdel=", "dr=", "eps=", "progdir=",
          "help" ] )
  except getopt.GetoptError, err:
    help()

  for o, a in opts:
    if o in ("--rmin",):
      rmin = float(a)
    elif o in ("--rmax",):
      rmax = float(a)
    elif o in ("--rdel", "--dr"):
      rdel = float(a)
    elif o in ("--eps",):
      epsv = float(a)
    elif o in ("--progdir",):
      progdir = a
    elif o in ("-h", "--help"):
      help()

  if len(args) > 0:
    fncfg = args[0]
  else:
    ls = glob.glob("*.cfg")
    if len(ls) == 0:
      print "please specify a .cfg file!"
      exit(1)
    else:
      fncfg = ls[0]

  if progdir == None:
    # search the path that contains the program
    progdir = os.path.join("rism", "prog")
    for i in range(10):
      if os.path.exists(progdir): break
      progdir = os.path.join("..", progdir)
    print "program directory: " + progdir

  prog = os.path.join(progdir, prog0)



def loadcfg(fn):
  ''' load the configuration template '''
  s = open(fn).readlines()
  n = len(s)
  # search for the last line that specifies the distance constraint
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
  fndump = "dump_r%s.dat" % r
  open(fncfg, "w").writelines(cfg)

  cmd = "%s %s %s > %s" % (prog, opts, fncfg, fndump)

  # run the command
  ret = os.system(cmd)
  if ret != 0:
    print "command [%s] failed with %s" % (cmd, ret)
    exit(1)

  # analyze the output
  bmu = 0
  if r != None:
    s = open(fndump).readlines()
    n = len(s)
    for i in range(n - 1, -1, -1):
      ln = s[i]
      m = re.search("mol.*chem.*pot.*X beta(.*)", ln)
      if m:
        bmu = float( m.group(1).strip() )
        break

  if r == None:
    os.rename("out.dat", name + "_out.dat")
  # remove unused files
  os.remove(fndump)
  os.remove(fncfg)
  os.remove("crdnum.dat")
  return bmu



def chempot():
  ''' compute the chemical potential '''
  # build the program
  os.system("make -C " + progdir)

  # load the configuration file
  cfgtemplate, cfgln = loadcfg(fncfg)
  name = fncfg.split(".")[0]

  i = int( rmin / rdel );
  r = i * rdel
  ls = []
  # scan over the distance r
  while r < rmax + .5 * rdel:
    # change the distance to r and run the script
    # retrieve the chemical potential
    bmu = runscript(cfgtemplate, cfgln, r)
    ls += [(r, bmu), ]
    i += 1
    r = i * rdel

  # run the script with infinite separation
  # this determines the base-line value of the chemical potential
  bmu0 = runscript(cfgtemplate, cfgln, 1e30)
  sout = "#    r    beta*mu\n"
  for i in range(len(ls)):
    sout += "%9s %s\n" % (ls[i][0], ls[i][1] - bmu0)

  # run the script without the distance constraint
  # this computes the chemical potential from t(r)
  runscript(cfgtemplate, cfgln, None, name)

  # save the PMF from the chemical potential route
  fnbmu = fncfg.split(".")[0] + "_bmu.dat"
  open(fnbmu, "w").write(sout)

  # prepare a gnuplot script
  fngp = name + "_pmf.gp"
  open(fngp, "w").write(r'''#!/usr/bin/env gnuplot
set encoding iso_8859_1
set terminal push
set terminal postscript eps enhanced size 5, 7 font "Times, 20"
set output "$NAME.eps"
set multiplot

# empirical dielectric constant
eps = $EPSV

unset xlabel
set ylabel "{/Symbol-Oblique b} {/Times-Italic W}^{ex}" offset 1, 0

set xrange [$RMIN:$RMAX]
set yrange [:]

set key spacing 1.5

ht = 0.5
hb = 1 - ht

set size 1, ht
set origin 0, hb

plot [:][:] \
  "$NAME_out.dat"  u 1:(($7 == 3 && $8 == 4)?-$3:1/0)  w l  lt 1       t "{/Symbol-Oblique b }{/Times-Italic W}^{ex}", \
  "$NAME_bmu.dat"  u 1:($2)                            w lp lt 4 pt 2  t "{/Symbol-Oblique b D m}^{ex}", \
  0 lt 1 lw 0.5 notitle

set size 1, hb
set origin 0, 0

set xlabel "{/Times-Italic r} ({\305})"
set ylabel "{/Symbol-Oblique b D} {/Times-Italic W}" offset 1, 0

plot [:][:] "$NAME_out.dat"  u 1:(($7 == 3 && $8 == 4)?$10+$9/eps:1/0) \
                             w l  lt 2       t "{/Symbol-Oblique b D}{/Times-Italic W}", \
  0 lt 1 lw 0.5 notitle

unset multiplot
unset output
set terminal pop
'''.replace("$NAME", str(name)
  ).replace("$RMIN", str(rmin)
  ).replace("$RMAX", str(rmax)
  ).replace("$EPSV", str(epsv)) )
  os.system("chmod 755 " + fngp)
  os.system("gnuplot " + fngp)



if __name__ == "__main__":
  doargs()
  chempot()
