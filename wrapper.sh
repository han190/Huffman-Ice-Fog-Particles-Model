#! /usr/bin/env bash

CURDIR=$(pwd)
BIN=$CURDIR/bin
OBJ=$CURDIR/obj
DAT=$CURDIR/data
SRC=$CURDIR/src

mkdir bin obj data

gfortran -o $OBJ/real.o -J$OBJ -c $SRC/real.f90
gfortran -o $OBJ/ice_particles.o -J$OBJ -c $SRC/ice_particles.f90
gfortran -o $OBJ/main.o -J$OBJ -c $SRC/main.f90
gfortran -o $BIN/ice_particles $OBJ/*.o

cd $DAT
../bin/ice_particles
gnuplot ../src/plot.plt
cd ..
