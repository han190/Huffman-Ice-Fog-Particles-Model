# A SIMPLE ICE FOG PARTICLES MODEL REBUILT FROM [HUFFMAN, 1957]

## Structure

This program contains five files, four fortran files and a gnuplot file:

1.  real.f90                    -- type and constant module
2.  ice_particles_eqn.f90       -- equation module
3.  ice_particles_solver.f90    -- solver module
4.  main.f90                    -- main program
5.  plot.plt                    -- plotting program

## Math

The main equation we are solving here is:

```tex
\frac{dS}{dt}
```

$dS/dt =  - L * M * S / (R * T^2) * dT/dt  - R * T / (e_s * M) * Sum(dm/dt)$

and the main idea to solve the equation is, for each time step i, assume the
time at step i is t_i, we have 

    L(t) => function L_i(t_i) 
    M    => constant
    S(t) => iteration S(t_i) = S(t_i-1) + dS/dt(t_i-1)
    R    => constant
    T(t) => function T_i(t_i)
    dT(t)/dt => function dT_i(t_i)/dt_i
    e_s(t) => es_i(t_i)
    r(t) => r_i(t_i) = r_i-1(t_i-1) + dr/dt(t_i-1)
    sdm(t)/dt => sdm(t_i)/dt = sdm(1 -> t_i-1)/dt