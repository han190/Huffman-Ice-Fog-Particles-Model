program main
    use real_m
    use ice_fog_particles_m
    use ice_fog_particles_solver_m
    implicit none 

    ! global variables
    real(dp) :: v_0, a, b, tt_i, tt_0, f
    type(ice_fog_particles_t), pointer :: ifp
    real(dp) :: t

    plot_temp: block 
        a = 5e-4_dp
        b = 66.7_dp
        v_0 = 2000._dp
        tt_i = 60._dp
        tt_0 = -40._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 1, file = "fig1_1.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.3e2_dp) exit
            write(1, *) t, ifp%tt(t) - tt_m
            t = t * 1.1_dp
        end do
        deallocate(ifp)

        a = 2e-4_dp
        b = 667._dp
        v_0 = 200._dp
        tt_i = 30._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 2, file = "fig1_2.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.3e2_dp) exit
            write(2, *) t, ifp%tt(t) - tt_m
            t = t * 1.1_dp
        end do
        deallocate(ifp)

        a = 5e-4_dp
        b = 6670._dp
        v_0 = 20._dp
        tt_i = 0._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 3, file = "fig1_3.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.3e2_dp) exit
            write(3, *) t, ifp%tt(t) - tt_m
            t = t * 1.1_dp
        end do
        deallocate(ifp)

    end block plot_temp

    plot_saturation_vapor_pres: block 
        a = 5e-4_dp
        b = 66.7_dp
        v_0 = 2000._dp
        tt_i = 80._dp
        tt_0 = -60._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 4, file = "fig2_1.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.e4_dp) exit
            write(4, *) ifp%tt(t), ifp%e_s(t)
            t = t * 1.1_dp
        end do
        deallocate(ifp)
    end block plot_saturation_vapor_pres

    plot_sigma_n_ll: block 
        a = 5e-4_dp
        b = 66.7_dp
        v_0 = 2000._dp
        tt_i = 80._dp
        tt_0 = -60._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 5, file = "fig3_1.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.e4_dp) exit 
            write (5, *) ifp%tt(t), ifp%ll(t)*1e-10_dp, ifp%sigma(t)
            t = t * 1.01_dp
        end do 
        deallocate(ifp)
    end block plot_sigma_n_ll

    plot_uv: block 
        a = 5e-4_dp
        b = 66.7_dp
        v_0 = 2000._dp
        tt_i = 80._dp
        tt_0 = -60._dp

        allocate(ifp)
        call ifp%eqn_init(v_0, a, b, tt_i, tt_0, f)
        open(unit = 6, file = "fig4_1.txt", status = "unknown")
        t = 1e-4_dp
        do
            if (t >= 1.e4_dp) exit 
            write (6, *) ifp%tt(t), ifp%u(t), ifp%v(t)
            t = t * 1.01_dp
        end do 
        deallocate(ifp)
    end block plot_uv

    ! solve_s: block 
    !     type(ice_fog_particles_solver_t), pointer :: ifps
    !     integer :: i 

    !     a = 5e-4_dp
    !     b = 66.7_dp
    !     v_0 = 2000._dp
    !     tt_i = 60._dp
    !     tt_0 = -40._dp
    !     f = 5e-4_dp

    !     allocate(ifps)
    !     call ifps%init(v_0, a, b, tt_i, tt_0, f)
    !     call ifps%slv()

    !     open(unit = 7, file = "fig7_1.txt", status = "unknown")
    !     print *, nmax
    !     print *, 'size = ', size(t_arr)
    !     write (7, *) "t", "r", "drdt", "S", "dSdt", "I"
    !     do i = 1, nmax
    !         write (7, *) t_arr(i), r_arr(i), drdt_arr(i), &
    !             ss_arr(i), dssdt_arr(i), ii_arr(i)
    !     end do 
    !     deallocate(ifps)
    ! end block solve_s
        
end program main 