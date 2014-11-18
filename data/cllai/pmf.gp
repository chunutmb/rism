#!/usr/bin/env gnuplot

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 20"
set output "NaCl.eps"

eps = 7

set xlabel "{/Times-Italic r} ({\305})"
set ylabel "{/Symbol-Oblique b} {/Times-Italic W}^{ex}" offset 1, 0

set xtics 2
set mxtics 2
set ytics 50
set mytics 5

set xrange [4:16]
set yrange [:]

set key spacing 1.5

plot [:][:] \
  "NaCl_out.dat"  u 1:(($7 == 3 && $8 == 4)?-$3:1/0)  w l  lt 1       t "{/Symbol-Oblique b }{/Times-Italic W}^{ex}", \
  "NaCl_bmu.dat"  u 1:($2)                            w lp lt 4 pt 2  t "{/Symbol-Oblique b D m}^{ex}", \
  0 lt 1 lw 0.5 notitle

unset output
set terminal pop

