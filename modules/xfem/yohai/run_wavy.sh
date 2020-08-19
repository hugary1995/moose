#!/usr/bin/env bash

for ((i = 0; i < 16; i++)); do
  echo ========================================================
  echo mpiexec -n 10 ../xfem-opt -i plate_wavy.i A=0.$(printf %02d $i)
  echo ========================================================
  mpiexec -n 10 ../xfem-opt -i plate_wavy.i A=0.$(printf %02d $i)
done
