#!/usr/bin/env gnuplot
#
# reproduction of Fig. 2 in
#
# Alkali halides in water: Ion-solvent correlations and ion-ion potentials of
# mean force at infinite dilution
# B. Montgomery Pettitt and Peter J. Rossky
# J. Chem. Phys. 84(10) 5836-5844 (1986)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 3.0 font "Times, 20"
set output "pr1986fig2.eps"

eps = 7

set title "Cl^- - H_2O"

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic g}({/Times-Italic r})" offset 1, 0

set xtics 2
set ytics 1

set key spacing 1.5

plot [0:8][0:5] \
  "outNaCl.dat" u 1:(($7 == 0 && $8 == 4)?1+$2+$3:1/0) w l lt 1 t "Cl^- - O", \
  "outNaCl.dat" u 1:(($7 == 1 && $8 == 4)?1+$2+$3:1/0) w l lt 2 t "Cl^- - H", \
  1 lt 4 lw 0.2 notitle

unset output
set terminal pop

