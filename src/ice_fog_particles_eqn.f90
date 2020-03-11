module ice_fog_particles_eqn_m
    use real_m
    use runge_kutta_m
    implicit none 
    private 

    type, public, extends(runge_kutta_t) :: ice_fog_particles_eqn_t
        real(real_t) :: v_0 ! initial velocity
        real(real_t) :: a, b ! coefficients required to compute temp_abs
        real(real_t) :: T_i ! initial temperature
        real(real_t) :: T_0 ! environment temperature
        real(real_t) :: E ! the decrease in water vapor density, usually 0
        real(real_t) :: M ! the molecular weight of water
        real(real_t) :: K ! the thermal conductivity of air 
        real(real_t) :: D ! the diffusivity of water vapor in air
        real(real_t) :: R ! the gas constant
        real(real_t) :: f ! kinetics of vapor molecules and particle surface
        real(real_t) :: e_i, e_w ! the saturation vapor pressure(ice and water)
        real(real_t) :: L_v, L_s ! the heats of vaporization and sublimation
    contains
        procedure :: abs_temp => abs_temp_func
        procedure :: derivs => derivs_sub
    end type ice_fog_particles_eqn_t

    real(real_t), parameter :: pi = 3.1415926535897932384626_real_t
    real(real_t), parameter :: T_m = 273._real_t ! K, the melting temperature
    real(real_t), parameter :: constA = 6.5e-4_real_t ! cm-3 sec-1, [Bigg, 1953]
    real(real_t), parameter :: eta = 1._real_t ! Avagrado's number 

contains

    subroutine init_sub( &
        this, v_0, a, b, &
        T_i, T_0, E, M, K, D, R, f, &
        e_i, e_w, L_v, L_s &
    )
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: v_0, a, b, T_i, T_0, E, M, K, D, R, f
        real(real_t), intent(in) :: e_i, e_w, L_v, L_s

        this%v_0 = v_0
        this%a, b = a, b
        this%T_i = T_i
        this%T_0 = T_0
        this%E = E
        this%M = M
        this%K = K
        this%D = D
        this%R = R
        this%f = f
        this%e_i = e_i
        this%e_w = e_w
        this%L_v = L_v
        this%L_s = L_s 
    end subroutine init_sub

    function abs_temp_func(this, t) result(temp)
        class(ice_fog_particles_eqn_t), intent(in) :: this
        real(real_t), intent(in) :: t ! time
        real(real_t) :: temp

        associate( &
            a => this%a, b => this%b, c => this%v_0/this%b, &
            T_i => this%T_i, T_0 => this%T_0 &
        )

            temp = T_0 + (T_i - T_0) / &
                (a * b * (T_i - T_0) * log(c * t + 1._real_t) + 1._real_t)
        end associate
    end function abs_temp_func 

    function dTdt_func(this, t) result(dTdt)
        class(ice_fog_particles_eqn_t), intent(in) :: this
        real(real_t), intent(in) :: t ! time 
        real(real_t) :: dTdt ! derivative of temperature

        associate( &
            a => this%a, b => this%b, c => this%v_0/this%b, &
            T_i => this%T_i, T_0 => this%T_0 &
        )
            dTdt = a * b * c * (T_i - T_0)**2 / (c * t + 1._real_t) / &
                (1._real_t + a * b * (T_i - T_0) * log(1._real_t + c * t))**2
        end associate
    end function dTdt_func

    function e_s_func(this, r, t) result(e_s)
        class(ice_fog_particles_eqn_t), intent(in) :: this 
        real(real_t), intent(in) :: r, t ! radius and time
        real(real_t) :: e_s, T_s

        T_s = T_m - abs_temp_func(this, t)
        associate(this%e_w => e_w, this%e_i => e_i)
            e_s = e_w + (e_i - e_w) * &
                ( 1._real_t - exp(-constA * r**3 * t * (exp(T_s) - 1._real_t)) )
        end associate
    end function e_s_func

    function L_func(this, r, t) result(L)
        class(ice_fog_particles_eqn_t), intent(in) :: this 
        real(real_t), intent(in) :: r, t ! radius and time
        real(real_t) :: L, T_s

        T_s = T_m - abs_temp_func(this, t)
        associate(this%L_v => L_v, this%L_s => L_s)
            L = L_v + (L_s - L_v) * &
                ( 1._real_t - exp(-constA * r**3 * t * (exp(T_s) - 1._real_t)) )
        end associate
    end function L_func

    function integral_term(this, t) result(I)
        class(ice_fog_particles_eqn_t), intent(in) :: this 
        real(real_t), intent(in) :: t 
        real(real_t) :: I 

    end function integral_term

    function sum_dmdt_func(this, S, r, t) result(dmdt)
        class(ice_fog_particles_eqn_t), intent(in) :: this
        real(real_t), intent(in) :: S, r, t
        real(real_t) :: dmdt 
        real(real_t) :: u, v 
        real(real_t) :: temp

        temp = abs_temp_func(this, t)

        associate( &
            L => this%L, M => this%M, K => this%K, R => this%R, &
            D => this%D, e_s => this%e_s &
        )
            u = L**2 * M / (K * R * temp**2)
            v = R * temp / (e_s * D * M)
            dmdt = 4._real_t * pi / (u + v) * (S - 1) * integral_term(this, t)
        end associate
    end function sum_dmdt_func

    subroutine derivs_sub(this, t, S, dSdt)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: t ! time
        real(real_t), dimension(:), intent(in) :: S ! the saturation ratio
        real(real_t), dimension(:), intent(out) :: dSdt
        real(real_t) :: temp, dTdt, e_s, dmdt

        dTdt = dTdt_func(this, t)
        temp = abs_temp_func(this, t)
        e_s = e_s_func(this, t)
        dmdt = sum_dmdt_func(this, s, r, t)
        associate( this%L => L, this%M => M, this%R => R, )
            dSdt(1) = - L*M*S / (R*temp**2) * dTdt - R*temp / (e_s*M)*dmdt
        end associate
    end subroutine derivs_sub

end module ice_fog_particles_eqn_m