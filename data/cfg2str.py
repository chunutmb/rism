#!/usr/bin/env python



''' convert a configuration file to a C-style string '''



import sys, os, re



def cfg2str(fn):
  ''' return the C string of the configuration file '''
  s = open(fn).readlines()
  n = len(s)
  for i in range(n):
    ln = s[i].strip()
    if ln == "" or ln.startswith("#"): s[i] = ""
    else: s[i] = re.sub(r"\s+", "", ln) + r"\n"
  return '"' + "".join(s) + '"'



if len(sys.argv) > 1:
  print cfg2str(sys.argv[1])
else:
  print "needs an configuration file"
