#!/usr/bin/env gnuplot
#
# reproduction of Fig. 8 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 5.0 font "Times, 20"
set output "pk1985fig8.eps"

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic g}({/Times-Italic r})" offset 1, 0

set xtics 2
set ytics 1

set xrange [0:8]

set key spacing 1.5

set style line 1 lt 1 lw 4
set style line 2 lt 2 lw 2
set style line 3 lt 4 lw 2
set style line 4 lt 5 lw 2
set style line 9 lt 4 lw 4 lc rgb "#cccccc"

set multiplot

ht = 0.5
hb = 1 - ht

set size 1, ht
set origin 0, hb

plot [:][0:] \
  "modelIIIc_out.dat" u 1:(($7 == 0 && $8 == 3)?1+$2+$3:1/0) w l ls 1 t "C^+ and O", \
  ""                  u 1:(($7 == 1 && $8 == 3)?1+$2+$3:1/0) w l ls 2 t "C^+ and H", \
  1 ls 9 notitle

set size 1, hb
set origin 0, 0

plot [:][0:] \
  "modelIIIc_out.dat" u 1:(($7 == 0 && $8 == 4)?1+$2+$3:1/0) w l ls 1 t "C^- and O", \
  ""                  u 1:(($7 == 1 && $8 == 4)?1+$2+$3:1/0) w l ls 2 t "C^- and H", \
  1 ls 9 notitle

unset multiplot

unset output
set terminal pop

