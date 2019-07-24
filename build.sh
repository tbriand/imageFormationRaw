#!/bin/bash
# Script for the compilation of the binaries

# binaries are put in the directory build/
dir=build/
mkdir $dir

# Compilation of the internal binaries
cd $dir
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

# Compilation of the external binaries
cd external/

# dcraw
gcc -o dcraw -O4 dcraw.c -lm -DNODEPS
mv dcraw ../$dir

# modified inverse compositional algorithm
cd modified_inverse_compositional
make
mv inverse_compositional_algorithm ../../$dir
make clean
cd ..

# ponomarenko
cd ponomarenko_v4/ponomarenko
make
mv ponomarenko ../../../$dir
make clean
cd ../../

cd ..

# move GNU parallel
cp external/parallel $dir
