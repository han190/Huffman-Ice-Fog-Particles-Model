module ice_particle_m
    use real_m
    implicit none
    private

    type, public :: ice_particle_t
        real(real_t) :: v_0, a, b
        real(real_t) :: T_i, T_0
        real(real_t) :: f
    contains
        procedure :: eqn_init => eqn_init_sub
        procedure :: temperature => temp_func
        procedure :: e_s => e_s_func
        procedure :: L => L_func
        procedure :: sigma => sigma_func
        procedure :: u => u_func, v => v_func
    end type ice_particle_t

    real(real_t), public, parameter :: pi = 3.1415926535897932384626_real_t
    real(real_t), public, parameter :: T_m = 273.15_real_t
    real(real_t), public, parameter :: constA = 6.5e-4_real_t
    real(real_t), public, parameter :: eta = 6.022e23_real_t
    real(real_t), public, parameter :: R = 8.31446261815324e7_real_t
    real(real_t), public, parameter :: M = 18.01528_real_t
    real(real_t), public, parameter :: K = 2.09e3_real_t
    real(real_t), public, parameter :: D0 = 0.106144_real_t

contains

    subroutine eqn_init_sub(this, v_0, a, b, T_i, T_0, f)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: v_0, a, b, T_i, T_0, f

        this%v_0 = v_0
        this%a   = a;   this%b   = b
        this%T_i = T_i; this%T_0 = T_0
        this%f   = f
    end subroutine eqn_init_sub

    function temp_func(this, t) result(temp)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: temp, deltaT, c

        associate(a=>this%a, b=>this%b, T_i=>this%T_i, T_0=>this%T_0)
            c = this%v_0/this%b
            deltaT = T_i - T_0
            temp = T_0 + deltaT / &
                (a*b*deltaT*log(c*t + 1._real_t) + 1._real_t)
        end associate
        temp = temp + T_m
    end function temp_func

    function dtemp_func(this, t) result(dTdt)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: dTdt, deltaT, c

        associate(a=>this%a, b=>this%b, T_i=>this%T_i, T_0=>this%T_0)
            c = this%v_0/this%b
            deltaT = T_i - T_0
            dTdt = - a*b*c*deltaT**2 / (1._real_t + c*t) / &
                (1._real_t + a*b*deltaT*log(1._real_t+c*t))**2
        end associate
    end function dtemp_func

    function e_w_func(this, t) result(e_w)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_w, temp
        real(real_t), parameter :: &
            C_1 = 8.72679635_real_t, &
            C_2 = 27.84317001_real_t, &
            C_3 = -17.82965432_real_t

        temp = temp_func(this, t)
        e_w = exp(C_1 + C_2*(1._real_t - T_m/temp) + C_3*log10(temp/T_m))
    end function e_w_func

    function e_i_func(this, t) result(e_i)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_i, temp
        real(real_t), parameter :: &
            C_1 = 8.91902223_real_t, &
            C_2 = -6.973508_real_t, &
            C_3 = -17.32425415_real_t

        temp = temp_func(this, t)
        e_i = exp(C_1 + C_2*(T_m/temp - 1._real_t) + C_3*log(T_m/temp))
    end function e_i_func

    function e_s_func(this, t) result(e_s)
        class(ice_particle_t), intent(in) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: e_s, e_w, e_i, temp

        temp = temp_func(this, t)
        e_w = e_w_func(this, t)
        e_i = e_i_func(this, t)

        if (temp > T_m) then
            e_s = e_w
        else
            e_s = e_w + (T_m - temp)/40._real_t*(e_i - e_w)
        end if
    end function e_s_func

    function L_func(this, t) result(L)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: L, temp

        temp = temp_func(this, t)

        if (temp > T_m) then
            L = -0.1069_real_t*(temp - T_m) + 71.5_real_t
        else
            L = -0.3129_real_t*(temp - 230._real_t) + 85._real_t
        end if

        L = ((L - 65._real_t)/40._real_t + 2.325_real_t)*1e10_real_t
    end function L_func

    function sigma_func(this, t) result(sigma)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: sigma, temp

        temp = temp_func(this, t)
        sigma = -0.1575_real_t*(temp - 220._real_t) + 83.9_real_t
    end function sigma_func

    function D_func(this, t) result(d) ! Emperical
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: d, temp

        temp = temp_func(this, t)
        d = D0*exp(0.007_real_t*(temp - 230._real_t))
    end function D_func

    function v_func(this, t) result(v)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: v, temp, e_s

        temp = temp_func(this, t)
        e_s = e_s_func(this, t)
        v = R * temp / (D_func(this, t) * M * e_s)
    end function v_func

    function u_func(this, t) result(u)
        class(ice_particle_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: u, temp, L

        temp = temp_func(this, t)
        L = L_func(this, t)
        u = L**2 * M / (K * R * temp**2)
    end function u_func

end module ice_particle_m