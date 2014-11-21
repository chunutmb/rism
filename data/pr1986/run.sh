#!/usr/bin/env bash

make -C ../../prog
prefix=pr1986
for tag in LiCl NaF NaCl KF KCl
do
  ../../prog/rism0 "$prefix""$tag".cfg -o out"$tag".dat
done
gnuplot "$prefix"fig1.gp
gnuplot "$prefix"fig2.gp
gnuplot "$prefix"fig3.gp
gnuplot "$prefix"fig4.gp
gnuplot "$prefix"fig5.gp


# for Fig. 6
python diffW.py "$prefix"NaCl
python diffW.py "$prefix"NaClz
gnuplot "$prefix"fig6.gp
