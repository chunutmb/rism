#!/usr/bin/env gnuplot
#
# reproduction of Fig. 1 in
#
# The Interionic potential of mean force in a molecular polar
# solvent from an extended RISM equation
# Fumio Hirata, Peter J. Rossky, and B. Montgomery Pettitt
# J. Chem. Phys. 78(6) 4133-4144 (1983)
#

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 8, 5.6 font "Times, 20"
set output "hrp1983fig1.eps"

set multiplot

wl = 0.55
wr = 1 - wl

ht = 0.48
hb = 1 - ht

set label 1 "(a)" at screen 0.09, 0.95 font "Times, 24"
set label 2 "(b)" at screen 0.09, 0.48 font "Times, 24"
set label 3 "(c)" at screen 0.57, 0.95 font "Times, 24"
set label 4 "(d)" at screen 0.57, 0.48 font "Times, 24"

set linestyle 1 lt 1 lw 4.0
set linestyle 2 lt 1 lw 2.0 lc rgb "#808080"
set linestyle 3 lt 2 lw 2.0 lc rgb "#404040"
set linestyle 4 lt 4 lw 0.5 lc rgb "#808080"

set xrange [2:8]
set yrange [0:30]

thexlabel = "{/Times-Italic r}({\305})"
theylabel = "{/Times-Italic N}({/Times-Italic r})"

set key Left reverse spacing 1.5

set size wl, ht
set origin 0, hb

unset xtics
set ylabel theylabel offset 1, 0
set ytics ("8." 8, "16." 16, "24." 24) nomirror
set lmargin 6
set rmargin 0
set tmargin 1
set bmargin 0

set key at 5.7, 29

plot [:][:] \
  "crdnumq0.5.dat"  u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 1 t "{/Times-Italic q} = 0.5{/Times-Italic e},  + +", \
  ""                u 1:(($3 == 2 && $4 == 1)?$2:1/0) w l ls 2 t "{/Times-Italic q} = 0.5{/Times-Italic e},  + -", \
  "crdnumq0.dat"    u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 3 t "{/Times-Italic q} = 0", \
  0 ls 4 notitle


set size wl, hb
set origin 0, 0

set xlabel thexlabel offset 0, 0.5
set xtics ("2.0" 2, "" 3, "4.0" 4, "" 5, "6.0" 6) nomirror
set format y "%g"
set tmargin 0
set bmargin 3

plot [:][:] \
  "crdnumq1.dat"    u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 1 t "{/Times-Italic q} = 1.0{/Times-Italic e},  + +", \
  ""                u 1:(($3 == 2 && $4 == 1)?$2:1/0) w l ls 2 t "{/Times-Italic q} = 1.0{/Times-Italic e},  + -", \
  "crdnumq0.dat"    u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 3 t "{/Times-Italic q} = 0", \
  0 ls 4 notitle



set size wr, ht
set origin wl, hb

set key at 6.0, 29

unset xtics
unset xlabel
unset ylabel
set ytics ("" 8, "" 16, "" 24)
set format y ""
set lmargin 0
set rmargin 1
set tmargin 1
set bmargin 0

plot [:][:] \
  "crdnumq2.dat"    u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 1 t "{/Times-Italic q} = 2.0{/Times-Italic e},  + +", \
  ""                u 1:(($3 == 2 && $4 == 1)?$2:1/0) w l ls 2 t "{/Times-Italic q} = 2.0{/Times-Italic e},  + -", \
  "crdnumq0.dat"    u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 3 t "{/Times-Italic q} = 0", \
  0 ls 4 notitle



set size wr, hb
set origin wl, 0

set key at 7.3, 29
set xlabel thexlabel offset 0, 0.5
set xtics ("2.0" 2, "" 3, "4.0" 4, "" 5, "6.0" 6) nomirror
unset ylabel
set ytics ("" 8, "" 16, "" 24)
set tmargin 0
set bmargin 3

plot [:][:] \
  "crdnumq1r.dat"   u 1:(($3 == 2 && $4 == 0)?$2:1/0) w l ls 1 t "{/Times-Italic q} = 1.0{/Times-Italic e}, {/Symbol-Oblique s} = 2.0{\305},  + +", \
  ""                u 1:(($3 == 2 && $4 == 1)?$2:1/0) w l ls 2 t "{/Times-Italic q} = 1.0{/Times-Italic e}, {/Symbol-Oblique s} = 2.0{\305},  + -", \
  0 ls 4 notitle



unset multiplot

unset output
set terminal pop

