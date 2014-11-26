#!/usr/bin/env gnuplot
#
# reproduction of Fig. 6 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 5.0 font "Times, 16"
set output "pk1985fig6.eps"

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic E} / {/Times-Italic kT}" offset 1, 0

set xtics 2.5
set ytics 1

set xrange [0:12]

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

plot [:][-1.3:1] \
  "modelIIIzfree_out.dat"   u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 1 t "Free Atoms", \
  "modelIIIzRe3.0_out.dat"  u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 2 t "{/Times-Italic R_e} = 3.0", \
  "modelIIIzRe2.5_out.dat"  u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 4 t "{/Times-Italic R_e} = 2.5", \
  "modelIIIzRe2.0_out.dat"  u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 3 t "{/Times-Italic R_e} = 2.0", \
  1 ls 9 notitle

set size 1, hb
set origin 0, 0

plot [:][-1.3:1] \
  "modelIIIzfree_out.dat"   u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 1 t "Free Atoms", \
  "modelIIIzRe1.40_out.dat" u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 2 t "{/Times-Italic R_e} = 1.40", \
  "modelIIIzRe1.54_out.dat" u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 4 t "{/Times-Italic R_e} = 1.54", \
  "modelIIIzRe1.75_out.dat" u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 3 t "{/Times-Italic R_e} = 1.75", \
  1 ls 9 notitle

unset multiplot

unset output
set terminal pop

