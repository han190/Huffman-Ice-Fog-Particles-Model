module real_m
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none 

    integer, parameter :: real_t = real64
    integer, public, parameter :: nmax = 10000 ! set the max
    real(real_t), public :: global_arr(5, nmax) ! [h, t, r, S, dSdt]

end module real_m 
