module real_m
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none 

    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    ! integer, public, parameter :: nmax = 10000
    ! real(dp), dimension(nmax), public :: t_arr, r_arr, drdt_arr
    ! real(dp), dimension(nmax), public :: ss_arr, dssdt_arr, ii_arr

end module real_m 
