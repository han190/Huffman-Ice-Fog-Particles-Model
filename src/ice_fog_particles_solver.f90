module ice_fog_particles_solver_m
    use real_m
    use ice_fog_particles_m
    implicit none 
    private 

    ! global arrays

    type, public :: ice_fog_particles_solver_t
        real(dp) :: v_0, a, b
        real(dp) :: tt_i, tt_0
        real(dp) :: f
        type(ice_fog_particles_t) :: ifp_ptr
    contains
        procedure :: init => solver_init_sub
        procedure :: slv => solve_sub
    end type ice_fog_particles_solver_t

contains

    subroutine solver_init_sub(this, v_0, a, b, tt_i, tt_0, f)
        class(ice_fog_particles_solver_t) :: this
        real(dp), intent(in) :: v_0, a, b, tt_i, tt_0, f

        !allocate(this%ifp_ptr)
        call this%ifp_ptr%eqn_init(v_0, a, b, tt_i, tt_0, f)
        
        this%v_0 = v_0
        this%a = a
        this%b = b
        this%tt_i = tt_i
        this%tt_0 = tt_0
        this%f = f

        r_arr = 0._dp
        drdt_arr = 0._dp
        ss_arr = 0._dp
        dssdt_arr = 0._dp
        ii_arr = 0._dp
    end subroutine solver_init_sub

    function integrand_func(this, step_i) result(integrand)
        class(ice_fog_particles_solver_t) :: this
        integer, intent(in) :: step_i
        real(dp) :: integrand, t, r, ii, ss

        t = t_arr(step_i)
        r = r_arr(step_i)
        ss = ss_arr(step_i)
        ii = this%ifp_ptr%ii(t, ss)
        integrand = ii * r**2 / (this%f + r)
    end function integrand_func

    function integral_func(this, step_i) result(integral)
        class(ice_fog_particles_solver_t) :: this
        integer, intent(in) :: step_i
        real(dp) :: integral
        integer :: i

        integral = 0._dp
        main_loop: do i = 1, step_i
            integral = integral + integrand_func(this, step_i)
        end do main_loop
    end function integral_func

    function sum_dmdt_func(this, step_i) result(sum_dmdt)
        class(ice_fog_particles_solver_t) :: this 
        integer, intent(in) :: step_i
        real(dp) :: sum_dmdt
        real(dp) :: integral, ss, u, v, t
        real(dp) :: top

        t = t_arr(step_i)
        u = this%ifp_ptr%u(t)
        v = this%ifp_ptr%v(t)
        ss = ss_arr(step_i)
        integral = integral_func(this, step_i)
        top = 4._dp * pi * (ss - 1._dp)
        sum_dmdt = top / (u + v) * integral
    end function sum_dmdt_func 

    function term_two_func(this, step_i) result(term_two)
        class(ice_fog_particles_solver_t) :: this 
        integer, intent(in) :: step_i
        real(dp) :: term_two
        real(dp) :: tt, e_s, t

        t = t_arr(step_i)
        tt = this%ifp_ptr%tt(t)
        e_s = this%ifp_ptr%e_s(t)
        term_two = rr * tt / (e_s * mm)
    end function term_two_func

    function term_one_func(this, step_i) result(term_one)
        class(ice_fog_particles_solver_t) :: this 
        integer, intent(in) :: step_i
        real(dp) :: term_one
        real(dp) :: dttdt, ss, tt, ll, t

        t = t_arr(step_i)
        tt = this%ifp_ptr%tt(t)
        ss = ss_arr(step_i)
        ll = this%ifp_ptr%ll(t)
        dttdt = this%ifp_ptr%dttdt(t)

        term_one = ll * mm * ss / (rr * tt**2) * dttdt
    end function term_one_func 

    function dssdt_func(this, step_i) result(dssdt)
        class(ice_fog_particles_solver_t) :: this 
        integer, intent(in) :: step_i
        real(dp) :: dssdt 
        real(dp) :: term_one, term_two, sum_dmdt

        term_one = term_one_func(this, step_i)
        term_two = term_two_func(this, step_i)
        sum_dmdt = sum_dmdt_func(this, step_i)

        dssdt =  - term_one - term_two * sum_dmdt
    end function dssdt_func

    subroutine solve_sub(this)
        class(ice_fog_particles_solver_t) :: this 
        integer :: i
        real(dp), parameter :: ss0 = 1._dp
        real(dp), parameter :: ii0 = 0._dp
        real(dp), parameter :: r0 = 0._dp
        real(dp), parameter :: t0 = 1.e-4_dp
        real(dp), parameter :: dt = 1.e-6_dp

        i = 1
        t_arr(1) = t0
        ss_arr(1) = ss0
        r_arr(1) = r0
        ii_arr(1) = ii0
        block
            real(dp) :: ll, a, b, c, deltatt, tt_i

            ll = this%ifp_ptr%ll(t_arr(1))
            a = this%a
            b = this%b 
            c = this%v_0 / this%b
            deltatt = this%tt_i - this%tt_0
            tt_i = this%tt_i
            dssdt_arr(1) = ll * mm * a * b * c * deltatt**2 / (rr * tt_i**2)
        end block
        drdt_arr(1) = 1.e-7_dp

        do
            if (i >= nmax) exit
            i = i + 1
            t_arr(i) = t_arr(i - 1) + dt
            ss_arr(i) = ss_arr(i - 1) + dssdt_arr(i - 1)
            ii_arr(i) = this%ifp_ptr%ii( t_arr(i), ss_arr(i) )
            r_arr(i) = r_arr(i - 1) + drdt_arr(i - 1)
            drdt_arr(i) = this%ifp_ptr%drdt( ss_arr(i), r_arr(i), t_arr(i) )
            dssdt_arr(i) = dssdt_func(this, i)
        end do
    end subroutine solve_sub

end module ice_fog_particles_solver_m
