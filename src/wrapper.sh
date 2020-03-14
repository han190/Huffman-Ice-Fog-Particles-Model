#! /usr/bin/env bash 

gfortran -c real.f90
gfortran -c rk_parameters.f90
gfortran -c runge_kutta.f90
gfortran -c ice_fog_particles_eqn.f90
gfortran -c ice_fog_particles_solver.f90
gfortran -c main.f90
gfortran -o ifp *.o
./ifp
gnuplot plot.plt
rm *.o *.mod *.txt
