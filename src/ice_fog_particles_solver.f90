module ice_fog_particles_solver_m
    use real_m
    use runge_kutta_m
    use ice_fog_particles_eqn_m
    implicit none
    private

    type, public :: ice_fog_particles_solver_t
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
        real(real_t) :: e_i, e_w ! the saturation vapor pressure(ice and water)
        real(real_t) :: L_v, L_s ! the heats of vaporization and sublimation
        type(ice_fog_particles_eqn_t), pointer :: ifp_eqn_p
    contains
        procedure :: init => ifp_solver_init_sub
    end type ice_fog_particles_solver_t

    real(real_t), parameter :: eps = 1.e-14_real_t
    real(real_t), parameter :: yscal = 1._real_t

contains

    subroutine ifp_solver_init_sub( &
        this, v_0, a, b, T_i, T_0, M, K, D, R, f, e_i, e_w, L_v, L_s &
    )
        class(ice_fog_particles_solver_t) :: this
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

        allocate(this%ifp_eqn_p)
        call this%ifp_eqn_p%rk_init(eps, yscal)
        call this%ifp_eqn_p%eqn_init( &
            v_0, a, b, T_i, T_0, M, K, &
            D, R, f, e_i, e_w, L_v, L_s &
        )
    end subroutine ifp_solver_init_sub

    subroutine solver_sub(this, tspan)
        class(ice_fog_particles_solver_t) :: this 
        real(real_t), intent(in) :: tspan(2)
        real(real_t) :: t0, tf, r0, h0
        real(real_t) :: tn, rn, Sn, hn, dSdtn
        real(real_t) :: hdid, hnext
        integer :: idx

        t0 = tspan(1); tf = tspan(2)
        r0 = 1.e-7_real_t
        h0 = 1.e-14_real_t
        global_arr = 0._real_t

        rn = r0; tn = t0; Sn = 1._real_t; hn = h0
        call this%ifp_eqn_p%derivs(tn, rn, Sn, dSdtn)
        global_arr(1:5, 1) = [hn, tn, rn, Sn, dSdtn]
        idx = 2

        do ! something 
            if (tn >= tf) exit 
            
            call this%ifp_eqn_p%derivs(tn, rn, Sn, dSdtn)
            call this%ifp_eqn_p%rkqs(rn, Sn, dSdtn, tn, hn, hdid, hnext)

            hn = hnext
            rn = rn + this%ifp_eqn_p%dr(Sn, rn, tn)
            idx = idx + 1
            global_arr(1:5, idx) = [hn, tn, rn, Sn, dSdtn]
        end do 
    end subroutine solver_sub

end module ice_fog_particles_solver_m
