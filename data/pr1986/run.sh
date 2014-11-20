#!/usr/bin/env bash

make -C ../../prog
prefix=pr1986
for tag in LiCl NaF NaCl KF KCl
do
  ../../prog/rism0 $prefix$tag.cfg -o $prefix$tag_out.dat > $prefix$tag_dump.dat
done
