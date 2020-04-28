program main
    use real_m
    use ice_fog_particles_solver_m
    implicit none 
    
    type(ice_fog_particles_solver_t), pointer :: ifps

    ! Curve A
    allocate(ifps)
    call ifps%solver_init( &
        2000._real_t, 5.e-4_real_t, &
        66.70_real_t, 60.00_real_t, &
        -44.0_real_t, 5.e-4_real_t &
    )
    print "(a)", "Compute saturation ratio for case A ..."
    call ifps%solver_slv( &
        "saturation_A.dat", &
        [1.e-4_real_t, 1.e-2_real_t] &
    )
    print "(a)", "Compute temperature for case A ..."
    call ifps%solver_plt_temp("temp_A.dat")
    deallocate(ifps)

    ! Curve B
    allocate(ifps)
    call ifps%solver_init( &
        200.0_real_t, 2.e-4_real_t, &
        667.0_real_t, 30.00_real_t, &
        -44.0_real_t, 5.e-4_real_t &
    )
    print "(a)", "Compute saturation ratio for case B ..."
    call ifps%solver_slv( &
        "saturation_B.dat", &
        [1.e-2_real_t, 1.e-0_real_t] &
    )
    print "(a)", "Compute temperature for case B ..."
    call ifps%solver_plt_temp("temp_B.dat")
    deallocate(ifps)

    ! Curve C
    allocate(ifps)
    call ifps%solver_init( &
        20.00_real_t, 5.e-4_real_t, &
        6670._real_t, 00.00_real_t, &
        -44.0_real_t, 5.e-4_real_t &
    )
    print "(a)", "Compute saturation ratio for case C ..."
    call ifps%solver_slv( &
        "saturation_C.dat", &
        [1.e-1_real_t, 1.e+1_real_t] &
    )
    print "(a)", "Compute temperature for case C ..."
    call ifps%solver_plt_temp("temp_C.dat")
    deallocate(ifps)
    
end program main
