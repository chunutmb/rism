#!/usr/bin/env bash

make -C ../../prog
for tag in q0 q0.5 q1 q2 q1r
do
  ../../prog/rism0 hrp1983$tag.cfg -o out$tag.dat -# crdnum$tag.dat
  gnuplot hrp1983$tag.gp
done
gnuplot hrp1983cmp.gp
gnuplot hrp1983crdnum.gp
#rm -f *.dat
