#!/usr/bin/env gnuplot
#
# reproduction of Fig. 5 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 2.5 font "Times, 20"
set output "pk1985fig5.eps"

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic E} / {/Times-Italic kT}" offset 2, 0

ymax = 0.7

set xtics 3
set ytics ymax

set xrange [0:12]

set key spacing 1.5

set title "Model III, nonpolar C_2 in TIPS water, {/Times-Italic R_e} = 1.7"

set style line 1 lt 1 lw 4
set style line 2 lt 2 lw 2
set style line 3 lt 4 lw 2
set style line 4 lt 5 lw 2
set style line 9 lt 4 lw 4 lc rgb "#cccccc"

plot [:][-ymax:ymax] \
  "modelIII_out.dat"  u 1:(($7 == 3 && $8 == 3)?$6-$3:1/0) w l ls 1 t "PMF", \
  ""                  u 1:(($7 == 3 && $8 == 3)?$6-$9:1/0) w l ls 2 t "{/Times-Italic U}^*", \
  1 ls 9 notitle

unset output
set terminal pop

