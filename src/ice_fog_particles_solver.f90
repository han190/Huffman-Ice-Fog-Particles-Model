module ice_fog_particles_solver_m
    use real_m
    use ice_fog_particles_eqn_m
    implicit none 
    private 

    type, public :: ice_fog_particles_solver_t
        real(dp) :: v_0, a, b
        real(dp) :: tt_i, tt_0
        real(dp) :: f
        type(ice_fog_particles_t), pointer :: ifp_ptr
    contains
        procedure :: init => solver_init_sub
        procedure :: slv => solve_sub
    end type ice_fog_particles_solver_t

contains

    subroutine solver_init_sub(this, v_0, a, b, tt_i, tt_0, f)
        class(ice_fog_particles_solver_t) :: this
        real(dp), intent(in) :: v_0, a, b, tt_i, tt_0, f

        allocate(this%ifp_ptr)
        call this%ifp_ptr%eqn_init(v_0, a, b, tt_i, tt_0, f)
        
        this%v_0 = v_0
        this%a = a
        this%b = b
        this%tt_i = tt_i
        this%tt_0 = tt_0
        this%f = f
    end subroutine solver_init_sub

    function dt_func(t) result(dt)
        real(dp), intent(in) :: t
        real(dp) :: dt

        if (t >= 1.e-6_dp .and. t < 1.e-3_dp) then
            dt = 1.e-5_dp
        else if (t >= 1.e-3_dp .and. t < 1.e-2_dp) then
            dt = 1.e-4_dp
        else if (t >= 1.e-2_dp .and. t < 1.e-1_dp) then
            dt = 1.e-3_dp
        else if (t >= 1.e-1_dp .and. t < 1._dp) then
            dt = 1.e-2_dp
        else if (t >= 1._dp .and. t < 10._dp) then
            dt = 1.e-1_dp
        else
            dt = 1._dp
        end if
    end function dt_func

    subroutine solve_sub(this, time, saturation_ratio)
        class(ice_fog_particles_solver_t) :: this
        real(dp), allocatable, intent(out) :: time(:), saturation_ratio(:)

        real(dp), parameter :: ss0 = 1._dp
        real(dp), parameter :: ii0 = 0._dp
        real(dp), parameter :: r0 = 0._dp
        real(dp), parameter :: t0 = 1.e-4_dp

        integer, parameter :: nmax = 100000
        real(dp), dimension(nmax) :: t_arr, ss_arr, integrand_arr

        real(dp) :: ss, ii, ll, tt, dttdt, e_s, u, v, r, t, dssdt, dt
        real(dp) :: integrand, integral
        real(dp) :: term_one, term_two, sum_dmdt
        real(dp) :: ss_min
        integer :: i

        ! Initialization
        dt              = 1.e-5_dp
        t_arr           = 0._dp
        ss_arr          = 0._dp
        integrand_arr   = 0._dp

        t   = t0
        ss  = ss0
        r   = r0
        ii  = ii0

        i = 1
        ss_min = 10._dp
        do
            if (i >= nmax .or. t >= 9.9) exit

            ! step i
            if (i /= 1) then
                r = this%ifp_ptr%r(t, ss)
                ii = this%ifp_ptr%ii(t, ss)
            end if
            tt          = this%ifp_ptr%tt(t)
            dttdt       = this%ifp_ptr%dttdt(t)
            u           = this%ifp_ptr%u(t)
            v           = this%ifp_ptr%v(t)
            e_s         = this%ifp_ptr%e_s(t)
            ll          = this%ifp_ptr%ll(t)
            integrand   = ii * r**2 / (this%f + r)
            integral    = sum(integrand_arr(1:i))

            term_one    = ll * mm * ss / (rr * tt**2) * dttdt
            term_two    = rr * tt / (e_s * mm)
            sum_dmdt    = 4._dp * pi * (ss - 1._dp) / (u + v) * integral
            dssdt       = - term_one - term_two * sum_dmdt

            if (dssdt < 0._dp) then
                if (ss_min >= ss) then
                    ss_min = ss
                else
                    exit
                end if
            end if

            ! store data to arrays
            t_arr(i)            = t
            ss_arr(i)           = ss
            integrand_arr(i)    = integrand

            ! step i+1
            i   = i + 1
            ss  = ss + dssdt*dt
            t   = t + dt
            dt  = dt_func(t)
        end do

        allocate(time(i), saturation_ratio(i))
        time(:) = t_arr(1:i)
        saturation_ratio(:) = ss_arr(1:i)
    end subroutine solve_sub

end module ice_fog_particles_solver_m
