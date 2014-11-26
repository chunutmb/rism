#!/usr/bin/env gnuplot
#
# reproduction of Fig. 3 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 5.0 font "Times, 20"
set output "pk1985fig3.eps"

set xtics 3
set ytics 0.5

set xrange [0:12]
set yrange [-1:0.5]

set key spacing 1.5

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 2
set style line 3 lt 4 lw 2
set style line 4 lt 5 lw 2
set style line 9 lt 4 lw 4 lc rgb "#cccccc"

set multiplot

ht = 0.47
hb = 1 - ht

set key bottom

dx = 0.27
dy = 0.08
set label 1 "(a)" at screen 1 - dx, 1 - dy - 0.04
set label 2 "(b)" at screen 1 - dx, hb - dy

set size 1, ht
set origin 0, hb
set lmargin 7
set bmargin 0

set format x ""
set ytics ("-0.5" -0.5, "0.0" 0.0, "0.5" 0.5)

plot [:][:] \
  "modelII_out.dat" u 1:(($7 == 2 && $8 == 3)?$6-$9:1/0)  w l ls 2 t "C_2, vacuum pair potential, {/Times-Italic U}^*({/Times-Italic r})", \
  "modelII_out.dat" u 1:(($7 == 2 && $8 == 3)?$6-$3:1/0)  w l ls 1 t "C_2, + -, PMF", \

set size 1, hb
set origin 0, 0

set tmargin 0
set bmargin 4

set format x "%g"
set xlabel "{/Times-Italic r}  in  {\305}"

set ytics ("-1.0" -1.0, "-0.5" -0.5, "0.0" 0.0, "0.5" 0.5)
set ylabel "{/Times-Italic E} / {/Times-Italic kT}" offset 2, 9

plot [:][:] \
  "modelII_out.dat" u 1:(($7 == 2 && $8 == 2)?$6-$9:1/0)  w l ls 2 t "C_2, vacuum pair potential, {/Times-Italic U}^*({/Times-Italic r})", \
  "modelII_out.dat" u 1:(($7 == 2 && $8 == 2)?$6-$3:1/0)  w l ls 1 t "C_2, + + or - -, PMF", \

unset multiplot

unset output
set terminal pop

