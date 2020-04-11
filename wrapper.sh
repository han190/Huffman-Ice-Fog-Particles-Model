#! /usr/bin/env bash

./wipper.sh 

CURDIR=$(pwd)
BIN=$CURDIR/bin
OBJ=$CURDIR/obj
DAT=$CURDIR/data
SRC=$CURDIR/src

mkdir bin obj data

gfortran -o $OBJ/real.o -J$OBJ -c $SRC/real.f90
gfortran -o $OBJ/ice_fog_particles_eqn.o -J$OBJ \
    -c $SRC/ice_fog_particles_eqn.f90
gfortran -o $OBJ/ice_fog_particles_solver.o -J$OBJ \
    -c $SRC/ice_fog_particles_solver.f90
gfortran -o $OBJ/main.o -J$OBJ -c $SRC/main.f90
gfortran -o $BIN/ice_fog_particles $OBJ/*.o

cd $DAT
../bin/ice_fog_particles
gnuplot ../src/plot.plt
# gnuplot ../src/plot4.plt
cd ..
