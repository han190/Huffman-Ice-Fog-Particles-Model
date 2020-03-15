module ice_fog_particles_eqn_m
    use real_m
    use runge_kutta_m
    implicit none 
    private 

    type, extends(runge_kutta_t), public :: ice_fog_particles_eqn_t
        real(real_t) :: v_0 ! initial velocity
        real(real_t) :: a, b ! coefficients required to compute temp_abs
        real(real_t) :: T_i ! initial temperature
        real(real_t) :: T_0 ! environment temperature
        !real(real_t) :: E ! the decrease in water vapor density, usually 0
        real(real_t) :: M ! the molecular weight of water
        real(real_t) :: K ! the thermal conductivity of air 
        real(real_t) :: D ! the diffusivity of water vapor in air
        real(real_t) :: R ! the gas constant
        real(real_t) :: f ! kinetics of vapor molecules and particle surface
        ! Don't if I should treat the following as an input or a function
        ! of time.
        real(real_t) :: e_i, e_w ! the saturation vapor pressure(ice and water)
        real(real_t) :: L_v, L_s ! the heats of vaporization and sublimation
    contains
        procedure :: eqn_init => eqn_init_sub
        procedure :: temp => temp_func
        procedure :: dr => drdt_func
        procedure :: r_embryo => I_func
        procedure :: derivs => derivs_sub
    end type ice_fog_particles_eqn_t

    real(real_t), public, parameter :: pi = 3.1415926535897932384626_real_t
    real(real_t), public, parameter :: T_m = 273.15_real_t
    ! K, the melting temperature
    real(real_t), public, parameter :: constA = 6.5e-4_real_t
    ! cm-3 sec-1, [Bigg, 1953]
    real(real_t), public, parameter :: eta = 6.022e23_real_t
    ! Avagrado's number mol^-1

contains

    subroutine eqn_init_sub( &
        this, v_0, a, b, T_i, T_0, M, K, D, R, f, e_i, e_w, L_v, L_s &
    )
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: v_0, a, b, T_i, T_0, M, K, D, R, f
        real(real_t), intent(in) :: e_i, e_w, L_v, L_s

        this%v_0 = v_0
        this%a   = a;   this%b   = b
        this%T_i = T_i; this%T_0 = T_0

        this%M   = M;   this%K   = K
        this%D   = D;   this%R   = R
        this%f   = f

        this%e_i = e_i; this%e_w = e_w
        this%L_v = L_v; this%L_s = L_s 
    end subroutine eqn_init_sub

    function temp_func(this, t) result(temp)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: t ! time
        real(real_t) :: temp, deltaT, c

        associate(a => this%a, b => this%b, T_i => this%T_i, T_0 => this%T_0)
            c = this%v_0/this%b
            deltaT = T_i - T_0
            temp = T_0 + deltaT / &
                (a*b*deltaT*log(c*t + 1._real_t) + 1._real_t)
        end associate
    end function temp_func

    function dtemp_func(this, t) result(dTdt)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: t
        real(real_t) :: dTdt, deltaT, c

        associate(a => this%a, b => this%b, T_i => this%T_i, T_0 => this%T_0)
            c = this%v_0/this%b
            deltaT = T_i - T_0
            dTdt = - a*b*c*deltaT**2 / (1._real_t + c*t) / &
                (1._real_t + a*b*deltaT*log(1._real_t+c*t))**2
        end associate
    end function dtemp_func

    function e_w_func(this, temp) result(e_w)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: temp
        real(real_t) :: e_w
        real(real_t), parameter :: e_s0 = 610.78_real_t ! Pascal
        real(real_t), parameter :: R_v = 461.5e4_real_t ! J/(kg K)
        real(real_t), parameter :: constL = 2.501e10 ! J/kg

        e_w = e_s0*exp(constL/R_v * (1._real_t/T_m - 1._real_t/temp))
        e_w = e_w*0.01_real_t
    end function e_w_func

    function e_i_func(this, temp) result(e_i)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: temp
        real(real_t) :: e_i
        real(real_t), parameter :: e_st = 611.20_real_t ! Pascal
        real(real_t), parameter :: R_v = 461.5e4_real_t ! J/(kg K)
        real(real_t), parameter :: constL = 2.834e10 ! J/kg

        e_i = e_st*exp(constL/R_v * (1._real_t/T_m - 1._real_t/temp))
        e_i = e_i*0.01_real_t
    end function e_i_func

    function e_s_func(this, r, t) result(e_s)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: r, t
        real(real_t) :: e_s, e_w, e_i, T_s, expt, temp

        temp = temp_func(this, t)
        ! e_w = e_w_func(this, temp)
        ! e_i = e_i_func(this, temp)

        if (temp >= 0._real_t) then
            T_s = 0._real_t
        else if (temp < 0._real_t) then
            T_s = temp
        end if

        expt = exp( - constA*r**3*t*(exp(T_s) - 1._real_t) )
        e_s  = e_w + (e_i - e_w) * ( 1._real_t - expt )
        associate(e_w => this%e_w, e_i => this%e_i)
            if (temp >= 0._real_t) then
                T_s = 0._real_t
            else if (temp < 0._real_t) then
                T_s = temp
            end if

            expt = exp( - constA*r**3*t*(exp(T_s) - 1._real_t) )
            e_s  = e_w + (e_i - e_w) * ( 1._real_t - expt )
        end associate
    end function e_s_func

    function L_func(this, r, t) result(lh)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: r, t 
        real(real_t) :: lh, T_s, expt, temp

        temp = temp_func(this, t)
        associate(L_s => this%L_s, L_v => this%L_v)
            if (temp >= 0._real_t) then
                T_s = 0._real_t
            else if (temp < 0._real_t) then
                T_s  = temp
            end if

            expt = exp( - constA*r**3*t*(exp(T_s) - 1._real_t) )
            lh   = L_s + (L_s - L_v) * (1._real_t - expt)
        end associate
    end function L_func

    subroutine UVcoeff_sub(this, r, t, u, v)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: r, t 
        real(real_t), intent(out) :: u, v
        real(real_t) :: temp, lh, e_s

        temp = temp_func(this, t)
        e_s  = e_s_func(this, r, t)
        lh   = L_func(this, r, t)

        u = lh**2 * this%M / (this%K * this%R * temp**2)
        v = this%R * temp / (e_s * this%D * this%M)
    end subroutine UVcoeff_sub

    function sigma_func(this, t) result(sigma)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: t 
        real(real_t) :: sigma, temp, delT
        real(real_t), parameter :: T_c = 647.15_real_t ! K
        real(real_t), parameter :: B = 235.8e-3_real_t ! N/m
        real(real_t), parameter :: lowerb = -0.625_real_t
        real(real_t), parameter :: mu = 1.256_real_t

        temp  = temp_func(this, t)
        delT  = (T_c - temp) / T_c
        sigma = B*delT**mu*(1._real_t + b*delT) * 1e-3
    end function sigma_func

    function drdt_func(this, S, r, t) result(drdt)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: S, r, t
        real(real_t) :: drdt, u, v
        real(real_t), parameter :: rho = .997_real_t

        call UVcoeff_sub(this, r, t, u, v)
        drdt = (S - 1._real_t) / ((u + v) * rho * r)
    end function drdt_func 

    function I_func(this, S, r, t) result(Rembryo)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: S, r, t
        real(real_t) :: Rembryo 
        real(real_t) :: sigma, temp, e_s
        real(real_t) :: term1, term2, expt

        sigma = sigma_func(this, t)
        temp  = temp_func(this, t)
        e_s   = e_s_func(this, r, t)

        associate(R => this%R, M => this%M)
            term1 = e_s**2 * S**2 / (R**2 * temp**2)
            term2 = sqrt(2._real_t * eta**3 * M * sigma / pi)
            expt  = exp( - 16._real_t * pi * eta * M**2 * sigma**3 / &
                (3._real_t * R**3 * temp**3 * log(S)**2 ) )
            Rembryo = term1 * term2 * expt
        end associate
    end function I_func

    subroutine integral_term_sub(this, S, r, t, ret)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: S, r, t 
        real(real_t), intent(out) :: ret
        real(real_t), parameter :: r_0 = 1.e-7_real_t
        real(real_t) :: Sn, rn, tn, tmp
        integer :: idx

        rn = r_0
        tmp = 0._real_t
        associate(f => this%f)
            do idx = 1, nmax
                tn = global_arr(2, idx)
                if (tn >= t .or. tn == 0._real_t) exit
                Sn = global_arr(3, idx)
                rn = global_arr(4, idx)

                tmp = tmp + I_func(this, Sn, rn, tn) * rn**2 / (f + rn)
            end do
        end associate

        ret = tmp
    end subroutine integral_term_sub

    function sdmdt_func(this, S, r, t) result(sdmdt)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: S, r, t
        real(real_t) :: sdmdt, u, v, intt

        call UVcoeff_sub(this, r, t, u, v)
        call integral_term_sub(this, S, r, t, intt)
        sdmdt = 4._real_t * pi / (u + v) * (S - 1) * intt
    end function sdmdt_func

    subroutine derivs_sub(this, t, r, S, dSdt)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: S, r, t
        real(real_t), intent(out) :: dSdt
        real(real_t) :: lh, e_s, sdmdt, temp, dTdt

        lh = L_func(this, r, t)
        e_s = e_s_func(this, r, t)
        sdmdt = sdmdt_func(this, S, r, t)
        temp = temp_func(this, t)
        dTdt = dtemp_func(this, t)

        associate(M => this%M, R => this%R)
            dSdt = - S/temp * (1._real_t - lh*M/(R*temp)) * dTdt &
                - R*temp/(e_s*M) * sdmdt
        end associate
    end subroutine derivs_sub

end module ice_fog_particles_eqn_m
