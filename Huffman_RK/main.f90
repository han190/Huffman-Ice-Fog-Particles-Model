program main
    use real_m
    use ice_fog_particles_solver_m
    implicit none 
    
    type(ice_fog_particles_solver_t), pointer :: ifps
    
    huffman_cases: block
        ! Curve A
        print "(a)", "compute curve A ..."; allocate(ifps)
        call ifps%solver_init(                                                 &
            2000._real_t, 5.e-4_real_t, 66.70_real_t,                          &
            60.00_real_t, -44.0_real_t, 5.e-4_real_t                           &
        )
        call ifps%solver_slv("fig1.txt", [1.e-4_real_t, 1.e-2_real_t])
        deallocate(ifps)
        ! Curve B
        print "(a)", "compute curve B ..."; allocate(ifps)
        call ifps%solver_init(                                                 &
            200.0_real_t, 2.e-4_real_t, 667.0_real_t,                          &
            30.00_real_t, -44.0_real_t, 5.e-4_real_t                           &
        )
        call ifps%solver_slv("fig2.txt", [1.e-2_real_t, 1.e-0_real_t])
        deallocate(ifps)
        ! Curve C
        print "(a)", "compute curve C ..."; allocate(ifps)
        call ifps%solver_init(                                                 &
            20.00_real_t, 5.e-4_real_t, 6670._real_t,                          &
            00.00_real_t, -44.0_real_t, 5.e-4_real_t                           &
        )
        call ifps%solver_slv("fig3.txt", [1.e-1_real_t, 1.e+1_real_t])
        deallocate(ifps)
    end block huffman_cases
    
end program main
