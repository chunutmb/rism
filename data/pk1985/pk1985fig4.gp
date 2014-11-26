#!/usr/bin/env gnuplot
#
# reproduction of Fig. 4 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 5.0 font "Times, 18"
set output "pk1985fig4.eps"

eps = 7

set xtics 3
set ytics 0.5

set xrange [0:12]
set yrange [0:1.0]

set key spacing 1.5

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 2
set style line 3 lt 4 lw 2
set style line 4 lt 5 lw 2
set style line 9 lt 4 lw 4 lc rgb "#cccccc"

set multiplot

ht = 0.33
hm = 0.30
hb = 1 - ht - hm

set key at 12, 0.8

dx = 0.27
dy = 0.04
set label 1 "(a)" at screen 1 - dx, 1 - dy - 0.03
set label 2 "(b)" at screen 1 - dx, hb + hm - dy
set label 3 "(c)" at screen 1 - dx, hb - dy

set size 1, ht
set origin 0, hb + hm
set lmargin 7
set bmargin 0

set format x ""
set ytics ("0.5" 0.5, "1.0" 1.0)

plot [:][:] \
  "modelIIasym_out.dat" u 1:(($7 == 2 && $8 == 2)?$6-$9+$11/eps:1/0)     w l ls 2 t "+ +, {/Times-Italic U}^* + {/Times-Italic U}'/{/Symbol-Oblique e}_{/Times-Italic p}", \
  "modelIIasym_out.dat" u 1:(($7 == 2 && $8 == 2)?$10+$6-$9+$11/eps:1/0) w l ls 1 t "+ +, PMF", \

set size 1, hm
set origin 0, hb

set tmargin 0
set bmargin 0

plot [:][:] \
  "modelIIasym_out.dat" u 1:(($7 == 3 && $8 == 3)?$6-$9+$11/eps:1/0)     w l ls 2 t "- -, {/Times-Italic U}^* + {/Times-Italic U}'/{/Symbol-Oblique e}_{/Times-Italic p}", \
  "modelIIasym_out.dat" u 1:(($7 == 3 && $8 == 3)?$10+$6-$9+$11/eps:1/0) w l ls 1 t "- -, PMF", \

set size 1, hb
set origin 0, 0

set bmargin 4

set format x "%g"
set xlabel "{/Times-Italic r}  in  {\305}"

set ytics ("0.0" 0.0, "0.5" 0.5, "1.0" 1.0)
set ylabel "{/Times-Italic E} / {/Times-Italic kT}" offset 1, 10

plot [:][:] \
  "modelIIasym_out.dat" u 1:(($7 == 2 && $8 == 3)?$6-$9+$11/eps:1/0)      w l ls 2 t "+ -, {/Times-Italic U}^* + {/Times-Italic U}'/{/Symbol-Oblique e}_{/Times-Italic p}", \
  "modelIIasym_out.dat" u 1:(($7 == 2 && $8 == 3)?$10+$6-$9+$11/eps:1/0)  w l ls 1 t "+ -, PMF", \

unset multiplot

unset output
set terminal pop

