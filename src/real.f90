module real_m
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none 

    integer, parameter :: sp = real64
    integer, parameter :: dp = real64
    integer, public, parameter :: nmax = 10000
    real(dp), dimension(nmax), public :: t_arr, r_arr, drdt_arr
    real(dp), dimension(nmax), public :: ss_arr, dssdt_arr, ii_arr
    ! integer, public, parameter :: nmax = 20000 ! set the max
    ! real(dp), public :: global_arr(5, nmax) ! [h, t, r, S, dSdt]

end module real_m 
