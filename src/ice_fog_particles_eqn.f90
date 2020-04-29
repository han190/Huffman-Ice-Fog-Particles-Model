module ice_fog_particles_eqn_m
    use real_m
    use runge_kutta_m
    implicit none
    private

    type, public, extends(runge_kutta_t) :: ice_fog_particles_t
        real(real_t) :: v_0, a, b, tt_i, tt_0, f
    contains
        procedure :: eqn_init => eqn_init_sub
        procedure :: tt => temp_func
        procedure :: dttdt => dtemp_func
        procedure :: e_s => e_s_func
        procedure :: ll => ll_func
        procedure :: sigma => sigma_func
        procedure :: u => u_func
        procedure :: v => v_func
        procedure :: ii => ii_func
        procedure :: r => r_func
        procedure :: derivs => derivs_sub
    end type ice_fog_particles_t

    real(real_t), public, parameter :: pi = 3.1415926535897932384626_real_t
    real(real_t), public, parameter :: tt_m = 273.15_real_t !check
    real(real_t), public, parameter :: eta = 6.022e23_real_t !check
    real(real_t), public, parameter :: rr = 8.31446261815324e7_real_t !check
    real(real_t), public, parameter :: mm = 18.01528_real_t !check
    real(real_t), public, parameter :: kk = 2.09e3_real_t

contains

    subroutine eqn_init_sub(this, v_0, a, b, tt_i, tt_0, f)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: v_0, a, b, tt_i, tt_0, f

        this%v_0 = v_0
        this%a   = a;   this%b   = b
        this%tt_i = tt_i; this%tt_0 = tt_0
        this%f   = f
    end subroutine eqn_init_sub

    function temp_func(this, t) result(temp)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: temp, deltaT, c

        associate(a=>this%a, b=>this%b, tt_i=>this%tt_i, tt_0=>this%tt_0)
            c = this%v_0/this%b
            deltaT = tt_i - tt_0
            temp = tt_0 + deltaT / &
                (a*b*deltaT*log(c*t + 1._real_t) + 1._real_t)
        end associate
        temp = temp + tt_m
    end function temp_func

    function dtemp_func(this, t) result(dTdt)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: dTdt, deltaT, c

        associate(a=>this%a, b=>this%b, tt_i=>this%tt_i, tt_0=>this%tt_0)
            c = this%v_0/this%b
            deltaT = tt_i - tt_0
            dTdt = - a*b*c*deltaT**2 / (1._real_t + c*t) / &
                (1._real_t + a*b*deltaT*log(1._real_t+c*t))**2
        end associate
    end function dtemp_func

    function e_w_func(this, t) result(e_w)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_w, temp
        real(real_t), parameter :: &
            c_1 = 8.72679635_real_t, &
            c_2 = 27.84317001_real_t, &
            c_3 = -17.82965432_real_t

        temp = temp_func(this, t)
        e_w = exp(c_1 + c_2*(1._real_t - tt_m/temp) + c_3*log10(temp/tt_m))
    end function e_w_func

    function e_i_func(this, t) result(e_i)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_i, temp
        real(real_t), parameter :: &
            c_1 = 8.9498907_real_t, &
            c_2 = -9.04755449_real_t, &
            c_3 = -36.00978668_real_t

        temp = temp_func(this, t)
        e_i = exp(c_1 + c_2*(tt_m/temp - 1._real_t) + c_3*log10(tt_m/temp))
    end function e_i_func

    function e_s_func(this, t) result(e_s)
        class(ice_fog_particles_t), intent(in) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_s, e_w, e_i, temp

        temp = temp_func(this, t)
        e_w = e_w_func(this, t)
        e_i = e_i_func(this, t)

        if (temp > tt_m) then
            e_s = e_w
        else
            e_s = e_w + (tt_m - temp)/40._real_t*(e_i - e_w)
        end if
    end function e_s_func

    function ll_func(this, t) result(ll)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: ll, temp

        temp = temp_func(this, t)

        if (temp > tt_m) then
            ll = -0.1069_real_t*(temp - tt_m) + 71.5_real_t
        else
            ll = -0.3129_real_t*(temp - 230._real_t) + 85._real_t
        end if

        ll = ((ll - 65._real_t)/40._real_t + 2.325_real_t)*1e10_real_t
    end function ll_func

    function sigma_func(this, t) result(sigma)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: sigma, temp

        temp = temp_func(this, t)
        sigma = -0.1575_real_t*(temp - 220._real_t) + 83.9_real_t
    end function sigma_func

    function dd_func(this, t) result(d) ! Emperical
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: d, temp
        real(real_t), parameter :: dd0 = 1.e-1_real_t
        real(real_t), parameter :: e_a = 5.e-1_real_t
        real(real_t) :: emperical_const

        temp = temp_func(this, t)
        emperical_const = (this%tt_0 + this%tt_i) + tt_m
        d = dd0*exp(-e_a*(temp - emperical_const))
    end function dd_func

    function v_func(this, t) result(v)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: v, temp, e_s

        temp = temp_func(this, t)
        e_s = e_s_func(this, t)
        v = rr * temp / (dd_func(this, t) * mm * e_s)
    end function v_func

    function u_func(this, t) result(u)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: u, temp, ll

        temp = temp_func(this, t)
        ll = ll_func(this, t)
        u = ll**2 * mm / (kk * rr * temp**2)
    end function u_func

    function ii_func(this, t, ss) result(ii)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t, ss
        real(real_t) :: ii, temp
        real(real_t) :: tt, bt ! top/bottom term
        real(real_t) :: expt ! exponential term
        real(real_t) :: term_1, term_2
        real(real_t) :: sigma, e_s

        temp = temp_func(this, t)
        sigma = sigma_func(this, t)
        e_s = e_s_func(this, t)

        tt = 16._real_t * pi * eta * mm**2 * sigma**3
        bt = 3._real_t * rr**3 * temp**3 * log(ss)**2
        expt = exp(- tt / bt)
        term_1 = (e_s * ss / (rr * temp))**2
        term_2 = sqrt(2._real_t * eta**3 * mm * sigma / pi)

        ii = term_1 * term_2 * expt
    end function ii_func

    function r_func(this, t, ss) result(r)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: t, ss
        real(real_t) :: r, r2, u, v

        u = u_func(this, t)
        v = v_func(this, t)
        r2 = 2._real_t * t * (ss - 1._real_t) / (u + v)
        
        if (r2 > 0._real_t) then 
            r = sqrt(r2)
        else
            r = 0._real_t
        end if
    end function r_func

    subroutine integrate(initial_condition, integrand, integral)
        logical, intent(in) :: initial_condition
        real(real_t), intent(in) :: integrand
        real(real_t), intent(out) :: integral 
        real(real_t), save :: temp

        if (initial_condition) then
            temp = integrand
        else
            temp = temp + integrand
        end if

        integral = temp
    end subroutine integrate

    subroutine derivs_sub(this, x, y, dydx)
        class(ice_fog_particles_t) :: this
        real(real_t), intent(in) :: x
        real(real_t), intent(in) :: y(:)
        real(real_t), intent(out) :: dydx(:)
        real(real_t) :: tt, dttdt, u, v, e_s, ll, ii, r
        real(real_t) :: integrand, integral
        real(real_t) :: term_one, term_two, sum_dmdt
        logical :: initial_condition

        initial_condition = x <= 1.e-5_real_t .and. y(1) == 1._real_t

        tt    = temp_func(this, x)
        dttdt = dtemp_func(this, x)
        u     = u_func(this, x)
        v     = v_func(this, x)
        e_s   = e_s_func(this, x)
        ll    = ll_func(this, x)
        ii    = ii_func(this, x, y(1))
        r     = r_func(this, x, y(1))

        integrand = ii * r**2 / (this%f + r)
        call integrate(initial_condition, integrand, integral)

        term_one = ll * mm * y(1) / (rr * tt**2) * dttdt
        term_two = rr * tt / (e_s * mm)
        sum_dmdt = 4._real_t * pi * (y(1) - 1._real_t) / (u + v) * integral
        dydx(1) = - term_one - term_two * sum_dmdt
    end subroutine derivs_sub

end module ice_fog_particles_eqn_m
