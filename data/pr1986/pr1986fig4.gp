#!/usr/bin/env gnuplot
#
# reproduction of Fig. 4 in
#
# Alkali halides in water: Ion-solvent correlations and ion-ion potentials of
# mean force at infinite dilution
# B. Montgomery Pettitt and Peter J. Rossky
# J. Chem. Phys. 84(10) 5836-5844 (1986)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 3.0 font "Times, 20"
set output "pr1986fig4.eps"

eps = 78

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic W}/{/Times-Italic kT}" offset 1, 0

set xtics 2
set ytics 1

set key reverse Left spacing 1.5

plot [2:8][0:4] \
  "outKCl.dat"  u 1:(($7 == 3 && $8 == 3)?$10+$6-$9+$9/eps:1/0) w l lt 1 t "K^+ K^+", \
  "outNaCl.dat" u 1:(($7 == 3 && $8 == 3)?$10+$6-$9+$9/eps:1/0) w l lt 2 t "Na^+ Na^+", \
  "outLiCl.dat" u 1:(($7 == 3 && $8 == 3)?$10+$6-$9+$9/eps:1/0) w l lt 5 t "Li^+ Li^+", \
  0 lt 4 lw 0.2 notitle

unset output
set terminal pop

