#!/usr/bin/env bash

make -C ../../prog
prefix=pk1985
prog=../../prog/rism0

for tag in modelI modelII modelIII
do
  $prog -V "$prefix""$tag".cfg -o "$tag"_out.dat
done

# change the densities to model I
rho1=0.00666666666666667
rho2=0.007
$prog -r1,"$rho1" -r2,"$rho1" -r3,"$rho1" "$prefix"modelI.cfg -o modelIrho0.02_out.dat
$prog -r1,"$rho2" -r2,"$rho2" -r3,"$rho2" "$prefix"modelI.cfg -o modelIrho0.021_out.dat


# change the partial charge to models II and III
$prog -q3,0.25 "$prefix"modelII.cfg -o modelIIasym_out.dat

$prog -q4,0.25 -q5,-0.25 "$prefix"modelIII.cfg -o modelIIIc_out.dat

# change the distance between the solute atoms
for tag in 2.0 2.5 3.0 1.75 1.54 1.40 Inf
do
  $prog -d4,5,"$tag" "$prefix"modelIII.cfg -o modelIIIRe"$tag"_out.dat
done

gnuplot "$prefix"fig1.gp
gnuplot "$prefix"fig2.gp
gnuplot "$prefix"fig3.gp
gnuplot "$prefix"fig4.gp
gnuplot "$prefix"fig5.gp
gnuplot "$prefix"fig6.gp
gnuplot "$prefix"fig7.gp
gnuplot "$prefix"fig8.gp


