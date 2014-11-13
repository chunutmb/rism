#!/usr/bin/env gnuplot
#
# reproduction of Fig. 2 in
#
# The Interionic potential of mean force in a molecular polar
# solvent from an extended RISM equation
# Fumio Hirata, Peter J. Rossky, and B. Montgomery Pettitt
# J. Chem. Phys. 78(6) 4133-4144 (1983)
#

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3, 5 font "Times, 20"
set output "hrp1983fig2.eps"

eps = 7

set xlabel "{/Times-Italic r} ({\305})"
set ylabel "{/Symbol-Oblique b} {/Times-Italic W}" offset 1, 0

set xtics 2,8,18
set mxtics 2
set ytics 0.5
set mytics 2

set xrange [2:18]
set yrange [-1:1]

set key spacing 1.5

# $6 is the total potential
# $9 is the electrostatic potential
# ($6 - $9) is the Lennard-Jones potential
# $10 is the short-range correction beta dW, see Eq. (50)

plot [:][:] \
  "outq0.dat" u 1:(($7 == 2 && $8 == 2)?$10+$6:1/0) w l lt 1 t "Total {/Symbol-Oblique b }{/Times-Italic W} (sum)", \
  "outq0.dat" u 1:(($7 == 2 && $8 == 2)?$6:1/0)     w l lt 4 t "Gas phase potential", \
  "outq0.dat" u 1:(($7 == 2 && $8 == 2)?$10:1/0)    w l lt 2 t "{/Symbol-Oblique b d}{/Times-Italic W}, Eq. (50)", \
  0 lt 1 lw 0.5 notitle




unset output
set terminal pop

