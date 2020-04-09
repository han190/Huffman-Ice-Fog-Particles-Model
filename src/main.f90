program main
    use real_m
    use ice_particle_m
    implicit none 

    ! Initializations
    real(real_t) :: v_0 = 0._real_t, a = 0._real_t, b = 0._real_t
    real(real_t) :: T_i = 0._real_t, T_0 = 0._real_t, f = 0._real_t

    plot_temp: block 
        type(ice_particle_t), pointer :: ice_particles
        real(real_t) :: t

        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 60._real_t
        T_0 = -40._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 1, file = "fig1_1.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.3e2_real_t) exit
            write(1, *) t, ice_particles%temperature(t) - T_m
            t = t * 1.1_real_t
        end do
        deallocate(ice_particles)

        a = 2e-4_real_t
        b = 667._real_t
        v_0 = 200._real_t
        T_i = 30._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 2, file = "fig1_2.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.3e2_real_t) exit
            write(2, *) t, ice_particles%temperature(t) - T_m
            t = t * 1.1_real_t
        end do
        deallocate(ice_particles)

        a = 5e-4_real_t
        b = 6670._real_t
        v_0 = 20._real_t
        T_i = 0._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 3, file = "fig1_3.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.3e2_real_t) exit
            write(3, *) t, ice_particles%temperature(t) - T_m
            t = t * 1.1_real_t
        end do
        deallocate(ice_particles)

    end block plot_temp

    plot_saturation_vapor_pres: block 
        type(ice_particle_t), pointer :: ice_particles
        real(real_t) :: t

        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 80._real_t
        T_0 = -60._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 4, file = "fig2_1.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.e4_real_t) exit
            write(4, *) ice_particles%temperature(t), &
                ice_particles%e_s(t)
            t = t * 1.1_real_t
        end do
        deallocate(ice_particles)
    end block plot_saturation_vapor_pres

    plot_sigma_n_L: block 
        type(ice_particle_t), pointer :: ice_particles
        real(real_t) :: t

        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 80._real_t
        T_0 = -60._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 5, file = "fig3_1.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.e4_real_t) exit 
            write (5, *) ice_particles%temperature(t), &
                ice_particles%L(t)*1e-10_real_t, ice_particles%sigma(t)
            t = t * 1.01_real_t
        end do 
        deallocate(ice_particles)
    end block plot_sigma_n_L

    plot_uv: block 
        type(ice_particle_t), pointer :: ice_particles
        real(real_t) :: t

        a = 5e-4_real_t
        b = 66.7_real_t
        v_0 = 2000._real_t
        T_i = 80._real_t
        T_0 = -60._real_t

        allocate(ice_particles)
        call ice_particles%eqn_init(v_0, a, b, T_i, T_0, f)
        open(unit = 6, file = "fig4_1.txt", status = "unknown")
        t = 1e-4_real_t
        do
            if (t >= 1.e4_real_t) exit 
            write (6, *) ice_particles%temperature(t), &
                ice_particles%u(t), ice_particles%v(t)
            t = t * 1.01_real_t
        end do 
        deallocate(ice_particles)
    end block plot_uv


        
end program main 