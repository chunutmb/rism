#!/usr/bin/env gnuplot
#
# reproduction of Fig. 6 in
#
# Alkali halides in water: Ion-solvent correlations and ion-ion potentials of
# mean force at infinite dilution
# B. Montgomery Pettitt and Peter J. Rossky
# J. Chem. Phys. 84(10) 5836-5844 (1986)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 6.0 font "Times, 20"
set output "pr1986fig6.eps"

set linestyle 1 lt 1 lw 4.0
set linestyle 2 lt 2 lw 2.0 lc rgb "#808080"
set linestyle 3 lt 5 lw 2.0 lc rgb "#404040"
set linestyle 4 lt 4 lw 0.5 lc rgb "#808080"

feps(t) = 87.740-0.40008*t+9.398e-4*t*t-1.410e-6*t*t*t
fdeps(t) = -0.4008+2*9.398e-4*t-3*1.410e-6*t*t

kB = 0.00198720414667
T = 298
T0 = 273.15
eps = feps(T-T0)
deps = fdeps(T-T0)
kT = kB*T

thexlabel = "{/Times-Italic r}  in  {\305}"
theylabel = "{/Times-Italic W} (kcal/mole)"

tW = "{/Times-Italic W}"
tS = "-{/Times-Italic T S}"
tE = "{/Times-Italic E}"

set xrange [2:8]

set xtics 2
set ytics 2

set key bottom reverse Left spacing 1.5

set multiplot

hb = 0.35
hm = 0.32
ht = 0.33

set label 1 "a)" at screen 0.92, hb-0.02
set label 2 "b)" at screen 0.92, hb+hm-0.02
set label 3 "c)" at screen 0.92, 1-0.04

set lmargin 7
set rmargin 1

set origin 0, 0
set size 1, hb

set tmargin 0
set xlabel thexlabel offset 0, 0.5
unset ylabel
set ytics ("-4.0" -4, "-2.0" -2, "0.0" 0, "2.0" 2)

plot [:][-4:4] \
  "pr1986NaCl_diff.dat"   u 1:2 w l ls 1 t tW, \
  ""                      u 1:3 w l ls 2 t tS, \
  ""                      u 1:4 w l ls 3 t tE, \
  0 ls 4 notitle



set origin 0, hb
set size 1, hm

set bmargin 0
unset xlabel
set format x ""
set ylabel theylabel offset 1, 0
set ytics ("-4.0" -4, "-2.0" -2, "0.0" 0, "2.0" 2)

plot [:][-4:4] \
  "outNaCl.dat"     u 1:(($7 == 3 && $8 == 4)?($6-$9+$9/eps)*kT:1/0)                    w l ls 1 t tW, \
  ""                u 1:(($7 == 3 && $8 == 4)?(-$9/eps/eps*T*deps*kT):1/0)              w l ls 2 t tS, \
  "outNaCl.dat"     u 1:(($7 == 3 && $8 == 4)?($6-$9+$9/eps+$9/eps/eps*T*deps)*kT:1/0)  w l ls 3 t tE, \
  0 ls 4 notitle



set origin 0, hm + hb
set size 1, ht

set tmargin 1
unset ylabel
set ytics ("-1.0" -1, "-0.5" -0.5, "0.0" 0, "0.5" 0.5)

plot [:][-1:0.5] \
  "pr1986NaClz_diff.dat"  u 1:2 w l ls 1 t tW, \
  ""                      u 1:3 w l ls 2 t tS, \
  ""                      u 1:4 w l ls 3 t tE, \
  0 ls 4 notitle



unset multiplot
unset output
set terminal pop

