#!/usr/bin/env bash

make -C ../../prog
prefix=pk1985
for tag in modelI modelIrho0.02 modelIrho0.021 modelII modelIIasym \
  modelIIIz modelIIIzRe2.0 modelIIIzRe2.5 modelIIIzRe3.0 \
  modelIIIzRe1.75 modelIIIzRe1.54 modelIIIzRe1.40 modelIIIzfree \
  modelIIIc
do
  ../../prog/rism0 -V "$prefix""$tag".cfg -o "$tag"_out.dat
done
gnuplot "$prefix"fig1.gp
gnuplot "$prefix"fig2.gp
gnuplot "$prefix"fig3.gp
gnuplot "$prefix"fig4.gp
gnuplot "$prefix"fig5.gp
gnuplot "$prefix"fig6.gp
gnuplot "$prefix"fig7.gp
gnuplot "$prefix"fig8.gp


