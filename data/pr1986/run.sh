#!/usr/bin/env bash

make -C ../../prog
prefix=pr1986
for tag in LiCl NaF NaCl KF KCl
do
  ../../prog/rism0 $prefix$tag.cfg -o out$tag.dat > dump$tag.dat
done
./"$prefix"fig1.gp
./"$prefix"fig2.gp
./"$prefix"fig3.gp
./"$prefix"fig4.gp
./"$prefix"fig5.gp
