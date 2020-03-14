program main
    use real_m; use runge_kutta_m
    use ice_fog_particles_eqn_m
    use ice_fog_particles_solver_m
    implicit none

    real(real_t) :: v_0 = 0._real_t
    real(real_t) :: a, b = 0._real_t
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

        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 60._real_t
        T_0 = -40._real_t

        M = 1._real_t !18.01528_real_t ! g/mol
        R = 1._real_t !8.314e7_real_t ! erg/(K mol)
        K = 1._real_t !0.02*1e5_real_t ! erg/(cm K)
        D = 1._real_t !0.242_real_t ! cm^2/s
        f = 1._real_t !5e-4_real_t ! 5 micro-meter
        e_i = 1._real_t !0.6_real_t
        e_w = 2._real_t !101._real_t
        L_v = 1._real_t ! just a random number
        L_s = 2._real_t ! just another random number

        allocate(ifp_slv_p)
        call ifp_slv_p%init(v_0,a,b,T_i,T_0,M,K,D,R,f,e_i,e_w,L_v,L_s)
        call ifp_slv_p%slv([1e-4_real_t, 1e2_real_t], "fig2_1.txt")
    end block plot_figure_2

end program main





