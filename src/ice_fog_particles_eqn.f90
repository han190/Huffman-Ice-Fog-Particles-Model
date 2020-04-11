module ice_fog_particles_m
    use real_m
    implicit none
    private

    ! C.G.ss. system is used.
    ! use double letter to represent upper case
    ! ignore division sign, for example: dr/dr => drdt

    type, public :: ice_fog_particles_t
        real(dp) :: v_0, a, b
        real(dp) :: tt_i, tt_0
        real(dp) :: f
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
        procedure :: drdt => drdt_func
    end type ice_fog_particles_t

    real(dp), public, parameter :: pi = 3.1415926535897932384626_dp !sure...
    real(dp), public, parameter :: tt_m = 273.15_dp !check
    real(dp), public, parameter :: eta = 6.022e23_dp !check
    real(dp), public, parameter :: rr = 8.31446261815324e7_dp !check
    real(dp), public, parameter :: mm = 18.01528_dp !check
    real(dp), public, parameter :: kk = 2.09e3_dp
    real(dp), public, parameter :: dd0 = 0.106144_dp

contains

    subroutine eqn_init_sub(this, v_0, a, b, tt_i, tt_0, f)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: v_0, a, b, tt_i, tt_0, f

        this%v_0 = v_0
        this%a   = a;   this%b   = b
        this%tt_i = tt_i; this%tt_0 = tt_0
        this%f   = f
    end subroutine eqn_init_sub

    function temp_func(this, t) result(temp)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: temp, deltaT, c

        associate(a=>this%a, b=>this%b, tt_i=>this%tt_i, tt_0=>this%tt_0)
            c = this%v_0/this%b
            deltaT = tt_i - tt_0
            temp = tt_0 + deltaT / &
                (a*b*deltaT*log(c*t + 1._dp) + 1._dp)
        end associate
        temp = temp + tt_m
    end function temp_func

    function dtemp_func(this, t) result(dTdt)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: dTdt, deltaT, c

        associate(a=>this%a, b=>this%b, tt_i=>this%tt_i, tt_0=>this%tt_0)
            c = this%v_0/this%b
            deltaT = tt_i - tt_0
            dTdt = - a*b*c*deltaT**2 / (1._dp + c*t) / &
                (1._dp + a*b*deltaT*log(1._dp+c*t))**2
        end associate
    end function dtemp_func

    function e_w_func(this, t) result(e_w)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: e_w, temp
        real(dp), parameter :: &
            c_1 = 8.72679635_dp, &
            c_2 = 27.84317001_dp, &
            c_3 = -17.82965432_dp

        temp = temp_func(this, t)
        e_w = exp(c_1 + c_2*(1._dp - tt_m/temp) + c_3*log10(temp/tt_m))
    end function e_w_func

    function e_i_func(this, t) result(e_i)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: e_i, temp
        real(dp), parameter :: &
            c_1 = 8.9498907_dp, &
            c_2 = -9.04755449_dp, &
            c_3 = -36.00978668_dp

        temp = temp_func(this, t)
        e_i = exp(c_1 + c_2*(tt_m/temp - 1._dp) + c_3*log10(tt_m/temp))
    end function e_i_func

    function e_s_func(this, t) result(e_s)
        class(ice_fog_particles_t), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: e_s, e_w, e_i, temp

        temp = temp_func(this, t)
        e_w = e_w_func(this, t)
        e_i = e_i_func(this, t)

        if (temp > tt_m) then
            e_s = e_w
        else
            e_s = e_w + (tt_m - temp)/40._dp*(e_i - e_w)
        end if
    end function e_s_func

    function ll_func(this, t) result(ll)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: ll, temp

        temp = temp_func(this, t)

        if (temp > tt_m) then
            ll = -0.1069_dp*(temp - tt_m) + 71.5_dp
        else
            ll = -0.3129_dp*(temp - 230._dp) + 85._dp
        end if

        ll = ((ll - 65._dp)/40._dp + 2.325_dp)*1e10_dp
    end function ll_func

    function sigma_func(this, t) result(sigma)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: sigma, temp

        temp = temp_func(this, t)
        sigma = -0.1575_dp*(temp - 220._dp) + 83.9_dp
    end function sigma_func

    function dd_func(this, t) result(d) ! Emperical
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: d, temp

        temp = temp_func(this, t)
        d = dd0*exp(0.007_dp*(temp - 230._dp))
    end function dd_func

    function v_func(this, t) result(v)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: v, temp, e_s

        temp = temp_func(this, t)
        e_s = e_s_func(this, t)
        v = rr * temp / (dd_func(this, t) * mm * e_s)
    end function v_func

    function u_func(this, t) result(u)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t
        real(dp) :: u, temp, ll

        temp = temp_func(this, t)
        ll = ll_func(this, t)
        u = ll**2 * mm / (kk * rr * temp**2)
    end function u_func

    function ii_func(this, t, ss) result(ii)
        class(ice_fog_particles_t) :: this
        real(dp), intent(in) :: t, ss
        real(dp) :: ii, temp
        real(dp) :: tt, bt ! top/bottom term
        real(dp) :: expt ! exponential term
        real(dp) :: term_1, term_2
        real(dp) :: sigma, e_s

        temp = temp_func(this, t)
        sigma = sigma_func(this, t)
        e_s = e_s_func(this, t)

        tt = 16._dp * pi * eta * mm**2 * sigma**3
        bt = 3._dp * rr**3 * temp**3 * log(ss)**2
        expt = exp(- tt / bt)
        term_1 = (e_s * ss / (rr * temp))**2
        term_2 = sqrt(2._dp * eta**3 * mm * sigma / pi)

        ii = term_1 * term_2 * expt
    end function ii_func

    function drdt_func(this, ss, r, t) result(drdt)
        class(ice_fog_particles_t) :: this 
        real(dp), intent(in) :: ss, r, t
        real(dp) :: drdt, u, v
        real(dp), parameter :: rho = .997_dp

        u = u_func(this, t)
        v = v_func(this, t)
        drdt = (ss - 1._dp) / ((u + v) * rho * r)
    end function drdt_func 

end module ice_fog_particles_m
