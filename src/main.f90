program main
    use real_m
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

        open(unit = 1, file="fig1_1.txt", status="unknown")
        t = 1e-4_real_t
        do
            if (t >= 1e2_real_t) exit 
            write(1, *) t, ifp_slv_p%ifp_eqn_p%temp(t)
            t = t * 1.1_real_t
        end do 
    end block plot_figure_1

end program main 





