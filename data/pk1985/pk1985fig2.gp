#!/usr/bin/env gnuplot
#
# reproduction of Fig. 2 in
#
# The potential of mean force between polyatomic molecules in polar
# molecular solvents
# B. Montgomery Pettitt and Martin Karplus
# J. Chem. Phys. 83(2) 781-789 (1985)

set encoding iso_8859_1

set terminal push
set terminal postscript eps enhanced size 3.5, 5.0 font "Times, 16"
set output "pk1985fig2.eps"

set xlabel "{/Times-Italic r}  in  {\305}"
set ylabel "{/Times-Italic E} / {/Times-Italic kT}" offset 1, 0

set xtics 2
set ytics 1

set xrange [0:12]
set yrange [-1:1]

set key spacing 1.5

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 2
set style line 3 lt 4 lw 2
set style line 4 lt 5 lw 2
set style line 9 lt 4 lw 4 lc rgb "#cccccc"

set multiplot

ht = 0.5
hb = 1 - ht

set size 1, ht
set origin 0, hb

plot [:][:] \
  "modelI_out.dat"          u 1:(($7 == 1 && $8 == 1)?$6:1/0)     w l ls 2 t "Ar_2-Ar_2 pair potential", \
  "modelI_out.dat"          u 1:(($7 == 1 && $8 == 1)?$6-$3:1/0)  w l ls 1 t "Ar_2-Ar_2 PMF, infinite dilution", \
  "modelIrho0.02_out.dat"   u 1:(($7 == 1 && $8 == 1)?$6-$3:1/0)  every 4 w p pt 7 ps 0.5 t "Ar_2-Ar_2 PMF, finite concentration"

set size 1, hb
set origin 0, 0

plot [:][:] \
  "modelIrho0.02_out.dat"   u 1:(($7 == 1 && $8 == 1)?$6-$3:1/0)  w l ls 2 t "Ar_2-Ar_2 PMF, {/Symbol-Oblqiue r} = 0.02", \
  "modelIrho0.021_out.dat"  u 1:(($7 == 1 && $8 == 1)?$6-$3:1/0)  w l ls 1 t "Ar_2-Ar_2 PMF, {/Symbol-Oblique r} = 0.021"

unset multiplot

unset output
set terminal pop

