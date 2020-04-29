module ice_fog_particles_solver_m
    use real_m
    use ice_fog_particles_eqn_m

    implicit none 
    private 

    type, public :: ice_fog_particles_solver_t
        real(real_t) :: v_0, a, b, tt_i, tt_0, f
        type(ice_fog_particles_t), pointer :: ice_fog_particles_p
    contains
        procedure :: solver_init => initialize_sub
        procedure :: solver_slv => solver_sub
        procedure :: solver_plt_temp => plot_temperature_sub
    end type ice_fog_particles_solver_t

    integer, parameter :: neq = 1
    integer, parameter :: nmax = 20000000 ! ten million points
    real(real_t), dimension(neq), parameter :: yscal = 1._real_t
    real(real_t), parameter :: eps = 1.e-16_real_t 

contains

    subroutine initialize_sub(this, v_0, a, b, tt_i, tt_0, f)
        class(ice_fog_particles_solver_t) :: this
        real(real_t), intent(in) :: v_0, a, b, tt_i, tt_0, f

        this%v_0 = v_0
        this%a = a
        this%b = b 
        this%tt_i = tt_i
        this%tt_0 = tt_0
        this%f = f

        allocate(this%ice_fog_particles_p)
        call this%ice_fog_particles_p%rk_init(eps, yscal)
        call this%ice_fog_particles_p%eqn_init(v_0, a, b, tt_i, tt_0, f)
    end subroutine initialize_sub

    subroutine solver_sub(this, filename, tspan)
        class(ice_fog_particles_solver_t) :: this
        character(len=*), intent(in) :: filename
        real(real_t), intent(in) :: tspan(2)
        real(real_t) :: ss(1), dssdt(1), t, htry, hdid, hnext
        real(real_t) :: ss_max, tol
        integer :: i 
        ! initial values
        real(real_t), parameter :: ss0 = 1._real_t, ii0 = 0._real_t
        real(real_t), parameter :: r0 = 0._real_t, t0 = 1.e-5_real_t

        t = t0
        ss = ss0
        htry = 1.e-14_real_t
        ss_max = 0._real_t
        
        open(unit = 1, file = filename, status = "unknown")
        
        i = 2
        main_loop: do
            if (i >= nmax .or. t >= tspan(2)) exit 

            call this%ice_fog_particles_p%derivs(t, ss, dssdt)
            call this%ice_fog_particles_p%rkqs(ss, dssdt, t, htry, hdid, hnext)
            htry = hnext
            
            if (t >= tspan(1)) then
                write (1, *) t, ss(1)
            end if 
            
            ss_max = max(ss(1), ss_max)
            if (ss_max > ss(1)) then
                tol = abs(ss(1) - 1._real_t)
                if (tol < t0) exit main_loop
            end if

            i = i + 1
        end do main_loop
    end subroutine solver_sub

    subroutine plot_temperature_sub(this, filename)
        class(ice_fog_particles_solver_t) :: this
        character(len=*), intent(in) :: filename
        real(real_t) :: t

        open(unit = 1, file = filename, status = "unknown")
        t = 1.e-4_real_t

        do
            if (t >= 1.1e2_real_t) exit 
            write (1, *) t, this%ice_fog_particles_p%tt(t) - tt_m
            t = t * 1.05_real_t
        end do
    end subroutine plot_temperature_sub


end module ice_fog_particles_solver_m

