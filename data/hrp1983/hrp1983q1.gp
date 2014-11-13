#!/usr/bin/env gnuplot
#
# reproduction of Fig. 4 in
#
# The Interionic potential of mean force in a molecular polar
# solvent from an extended RISM equation
# Fumio Hirata, Peter J. Rossky, and B. Montgomery Pettitt
# J. Chem. Phys. 78(6) 4133-4144 (1983)
#

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 5, 5 font "Times, 20"
set output "hrp1983fig4.eps"

eps = 7

set multiplot

wl = 0.55
wr = 1 - wl

set size wl, 1
set origin 0, 0
set title "+ +"

set xlabel "{/Times-Italic r} ({\305})"
set ylabel "{/Symbol-Oblique b} {/Times-Italic W}" offset 1, 0

set xtics 2,8,18
set mxtics 2
set ytics 45
set mytics 2

set xrange [2:18]
set yrange [-90:90]

set key spacing 1.5

# $6 is the total potential
# $9 is the electrostatic potential
# ($6 - $9) is the Lennard-Jones potential
# $10 is the short-range correction beta dW, see Eq. (50)

plot [:][:] \
  "outq1.dat" u 1:(($7 == 2 && $8 == 2)?$10+$6-$9+$9/eps:1/0) w l lt 1 t "{/Symbol-Oblique b }{/Times-Italic W}", \
  "outq1.dat" u 1:(($7 == 2 && $8 == 2)?$10:1/0)              w l lt 2 t "{/Symbol-Oblique b d}{/Times-Italic W}, Eq. (50)", \
  "outq1.dat" u 1:(($7 == 2 && $8 == 2)?$6-$9+$9/eps:1/0)     w l lt 4 t "{/Symbol-Oblique b }{/Times-Italic W_c}, Eq. (49)", \
  "outq0.dat" u 1:(($7 == 2 && $8 == 2)?$10:1/0) every 20     w p pt 2 t "{/Symbol-Oblique b d}{/Times-Italic W}^0, Eq. (50)", \
  0 lt 1 lw 0.5 notitle




set size wr, 1
set origin wl, 0
set title "+ -"

unset ylabel
set format y ""

plot [:][:] \
  "outq1.dat" u 1:(($7 == 2 && $8 == 3)?$10+$6-$9+$9/eps:1/0) w l lt 1 t "{/Symbol-Oblique b }{/Times-Italic W}", \
  "outq1.dat" u 1:(($7 == 2 && $8 == 3)?$10:1/0)              w l lt 2 t "{/Symbol-Oblique b d}{/Times-Italic W}, Eq. (50)", \
  "outq1.dat" u 1:(($7 == 2 && $8 == 3)?$6-$9+$9/eps:1/0)     w l lt 4 t "{/Symbol-Oblique b }{/Times-Italic W_c}, Eq. (49)", \
  "outq0.dat" u 1:(($7 == 2 && $8 == 2)?$10:1/0) every 20     w p pt 2 t "{/Symbol-Oblique b d}{/Times-Italic W}^0, Eq. (50)", \
  0 lt 1 lw 0.5 notitle

unset multiplot

unset output
set terminal pop

