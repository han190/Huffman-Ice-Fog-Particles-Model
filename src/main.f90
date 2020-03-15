program main
    use real_m; use runge_kutta_m
    use ice_fog_particles_eqn_m
    use ice_fog_particles_solver_m
    implicit none

    ! Initializations
    real(real_t) :: v_0 = 0._real_t
    real(real_t) :: a = 0._real_t, b = 0._real_t
    real(real_t) :: T_i = 0._real_t
    real(real_t) :: T_0 = 0._real_t
    real(real_t) :: M = 0._real_t
    real(real_t) :: K = 0._real_t
    real(real_t) :: D = 0._real_t
    real(real_t) :: R = 0._real_t
    real(real_t) :: f = 0._real_t
    real(real_t) :: e_i = 0._real_t, e_w = 0._real_t
    real(real_t) :: L_v = 0._real_t, L_s = 0._real_t

    plot_figure_1: block
        type(ice_fog_particles_solver_t), pointer :: ifp_slv_p
        real(real_t) :: t
        
        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 60._real_t
        T_0 = -40._real_t

        allocate(ifp_slv_p)
        call ifp_slv_p%init(v_0,a,b,T_i,T_0,M,K,D,R,f,e_i,e_w,L_v,L_s)

        open(unit = 1, file = "fig1_1.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1e2_real_t) exit
            write(1, *) t, ifp_slv_p%ifp_eqn_p%temp(t)
            t = t * 1.1_real_t
        end do
        deallocate(ifp_slv_p)

        a = 2e-4_real_t
        b = 667._real_t
        v_0 = 200._real_t
        T_i = 30._real_t

        allocate(ifp_slv_p)
        call ifp_slv_p%init(v_0,a,b,T_i,T_0,M,K,D,R,f,e_i,e_w,L_v,L_s)

        open(unit = 2, file = "fig1_2.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1e2_real_t) exit
            write(2, *) t, ifp_slv_p%ifp_eqn_p%temp(t)
            t = t * 1.1_real_t
        end do
        deallocate(ifp_slv_p)

        a = 5e-4_real_t
        b = 6670._real_t
        v_0 = 20._real_t
        T_i = 0._real_t

        allocate(ifp_slv_p)
        call ifp_slv_p%init(v_0,a,b,T_i,T_0,M,K,D,R,f,e_i,e_w,L_v,L_s)

        open(unit = 3, file = "fig1_3.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1e2_real_t) exit
            write(3, *) t, ifp_slv_p%ifp_eqn_p%temp(t)
            t = t * 1.1_real_t
        end do
        deallocate(ifp_slv_p)

    end block plot_figure_1

    plot_figure_2: block
        type(ice_fog_particles_solver_t), pointer :: ifp_slv_p
        real(real_t) :: tspan, arr(5, nmax)
        integer :: idx

        M   = 18._real_t !18._real_t ! g/mol
        R   = 8.314e7_real_t !8.314e7_real_t ! erg/(K mol)
        K   = 2.e3_real_t ! erg/(cm K)
        D   = .242_real_t ! cm^2/s
        f   = 5.e-4_real_t ! 5 micro-meter
        e_i = .6_real_t!0._real_t !.5_real_t ! 0.6_real_t
        e_w = 101._real_t!90._real_t !2._real_t ! 101._real_t
        L_v = 44000.e7_real_t
        L_s = 1e4_real_t

        a   = 5e-4_real_t
        b   = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 60._real_t
        T_0 = -40._real_t

        allocate(ifp_slv_p)
        call ifp_slv_p%init(v_0,a,b,T_i,T_0,M,K,D,R,f,e_i,e_w,L_v,L_s)
        call ifp_slv_p%slv([1e-4_real_t, 1e2_real_t], "fig2_1.txt")
    end block plot_figure_2

end program main
