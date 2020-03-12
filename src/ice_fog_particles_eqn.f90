module ice_fog_particles_eqn_m
    use real_m
    implicit none 
    private 

    type, public :: ice_fog_particles_eqn_t
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
        real(real_t) :: eta 
    contains
        procedure :: temperature => temp_func
    end type ice_fog_particles_eqn_t

    real(real_t), parameter :: pi = 3.1415926535897932384626_real_t
    real(real_t), parameter :: T_m = 273._real_t ! K, the melting temperature
    real(real_t), parameter :: constA = 6.5e-4_real_t ! cm-3 sec-1, [Bigg, 1953]
    real(real_t), parameter :: eta = 1._real_t ! Avagrado's number 

contains

    subroutine init_sub( &
        this, v_0, a, b, T_i, T_0, E, M, K, D, R, f, e_i, e_w, L_v, L_s &
    )
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: v_0, a, b, T_i, T_0, E, M, K, D, R, f
        real(real_t), intent(in) :: e_i, e_w, L_v, L_s

        this%v_0 = v_0
        this%a   = a;   this%b   = b
        this%T_i = T_i; this%T_0 = T_0
        this%E   = E

        this%M   = M;   this%K   = K
        this%D   = D;   this%R   = R
        this%f   = f

        this%e_i = e_i; this%e_w = e_w
        this%L_v = L_v; this%L_s = L_s 
    end subroutine init_sub

    function temp_func(this, t) result(temp)
        class(ice_fog_particles_eqn_t) :: this
        real(real_t), intent(in) :: t ! time
        real(real_t) :: temp, deltaT

        associate( &
            this%a => a, this%b => b, this%v_0/this%b => c, &
            this%T_i = T_i, this%T_0 = T_0 &     
        )
            
            deltaT = T_i - T_0
            temp = T_0 + deltaT / &
                (a*b*deltaT*log(c*t + 1._real_t) + 1._real_t)
        end associate
    end function temp_func

    function e_s_func(this, r, t) result(e_s)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: r, t
        real(real_t) :: e_s, T_s, expt, temp

        temp = temp_func(this, t)
        associate(this%e_w => e_w, this%e_i => e_i)
            T_s  = T_m - temp
            expt = exp( - constA*r**3*t*(exp(T_s) - 1._real_t) )
            e_s  = e_w + (e_i - e_w) * ( 1._real_t - expt )
        end associate
    end function e_s_func

    function L_func(this, r, t) result(lh)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: r, t 
        real(real_t) :: lh, T_s, expt, temp

        temp = temp_func(this, t)
        associate(this%L_s => L_s, this%L_v => L_v)
            T_s  = T_m - temp
            expt = exp( - constA*r**3*t*(exp(T_s) - 1._real_t) )
            lh   = L_s + (L_s - L_v) * ( 1._real_t - expt )
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
        v = this%R * this%T / (e_s * this%D * this%M)
    end subroutine UVcoeff_sub

    function sigma_func(this, t) result(sigma)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: t 
        real(real_t) :: sigma, temp, term1, term2
        real(real_t), parameter :: T_c = 647.15_real_t ! K
        real(real_t), parameter :: B = 235.8e-3_real_t ! N/m
        real(real_t), parameter :: lowerb = -0.625_real_t
        real(real_t), parameter :: mu = 1.256_real_t

        temp  = temp_func(this, t)
        term1 = ((T_c - temp) / T_c)**mu 
        term2 = (1._real_t + lowerb*((T_c - temp) / T_c))
        sigma = B*term1*term2
    end function sigma_func

    function drdt_func(this, S, r, t) result(drdt)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: S, r, t
        real(real_t) :: drdt, u, v
        real(real_t), parameter :: rho = 1._real_t

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

        associate( this%eta => eta, this%R => R, this%M => M )
            term1 = e_s**2 * S**2 / (R**2 * temp**2)
            term2 = sqrt(2._real_t * eta**3 * M * sigma / pi)
            expt  = exp( - 16._real_t * pi * eta * M**2 * sigma**3 / &
                (3._real_t * R**3 * temp**3 * log(S)**2 ) )
            Rembryo = term1 * term2 * expt
        end associate
    end function I_func

    subroutine integral_t(this, S, r, t, int_t)
        class(ice_fog_particles_eqn_t) :: this 
        real(real_t), intent(in) :: S, r, t 
        real(real_t), intent(out) :: int_t 
        real(real_t), parameter :: r_0 = 1.e-7_real_t
        real(real_t) :: r_i, r_ip1, t_i, t_ip1

        r_i = r_0
        associate(this%f => f)
            do
                r_ip1 = drdt(this, S, )
                
            end do
        end associate

end module ice_fog_particles_eqn_m