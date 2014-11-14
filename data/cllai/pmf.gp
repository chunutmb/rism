#!/usr/bin/env gnuplot
#

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 7, 5 font "Times, 20"
set output "neutral.eps"

eps = 7

set xlabel "{/Times-Italic r} ({\305})"
set ylabel "{/Symbol-Oblique b} {/Times-Italic W}^{ex}" offset 1, 0

set xtics 2
set mxtics 2
set ytics 0.1
set mytics 2

set xrange [4:16]
set yrange [:]

set key spacing 1.5

# $6 is the total potential
# $9 is the electrostatic potential
# ($6 - $9) is the Lennard-Jones potential
# $10 is the short-range correction beta dW, see Eq. (50)

plot [:][:] \
  "out0.dat"  u 1:(($7 == 3 && $8 == 4)?-$3:1/0)  w l  lt 1       t "{/Symbol-Oblique b }{/Times-Italic W}^{ex}", \
  "bmu.dat"   u 1:($2-28.4902)                    w lp lt 4 pt 2  t "{/Symbol-Oblique b D m}^{ex}", \
  0 lt 1 lw 0.5 notitle




unset output
set terminal pop

