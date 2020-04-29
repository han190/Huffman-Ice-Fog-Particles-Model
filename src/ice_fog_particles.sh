#! /usr/bin/env bash

echo "compile real.f90 ..."
gfortran -c real.f90
echo "compile rk_parameters.f90 ..."
gfortran -c rk_parameters.f90
echo "compile runge_kutta.f90 ..."
gfortran -c runge_kutta.f90
echo "compile ice_fog_particles_eqn.f90 ..."
gfortran -c ice_fog_particles_eqn.f90
echo "compile ice_fog_particles_solver.f90 ..."
gfortran -c ice_fog_particles_solver.f90
echo "compile main.f90 ..."
gfortran -c main.f90
echo "generate executables ..."
gfortran -o ice_fog_particles *.o
echo "execute binary files ..."
./ice_fog_particles
echo "plot output files ..."
gnuplot plot.plt
