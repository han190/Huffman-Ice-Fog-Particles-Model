module runge_kutta_m
    use real_m 
    implicit none 
    private

    integer, parameter :: nmax = 10000
    real(real_t) :: Sarr_global(4, nmax) ! [h, t, S, dSdt] 

    type, abstract, public :: runge_kutta_t 
        real(real_t) :: eps, yscal
    contains 
        procedure :: rk_init => initialize_sub 
        procedure :: rkck => rkck_sub 
        procedure :: rkqs => rkqs_sub 
        procedure(derivsx), deferred :: derivs 
    end type

    abstract interface
        subroutine derivsx(this, idx, t, S, Sarr, dSdt)
            use real_m
            import :: runge_kutta_t
            implicit none

            class(runge_kutta_t) :: this
            integer, intent(in) :: idx ! index
            real(real_t), intent(in) :: t, S
            real(real_t), intent(inout) :: Sarr(:,:)
            real(real_t), intent(out) :: dSdt
        end subroutine derivsx
    end interface

contains 

    subroutine initialize_sub(this, eps, yscal)
        class(runge_kutta_t) :: this 
        real(real_t), intent(in) :: eps, yscal

        this%eps   = eps 
        this%yscal = yscal 
    end subroutine initialize_sub 

    subroutine rkck_sub(this, idx, S, Sarr, dSdt, t, h, Sout, Serr)
        use rk_parameters_m
        implicit none 

        class(runge_kutta_t) :: this
        integer, intent(in) :: idx
        real(real_t), intent(in) :: t, S, dSdt, h
        real(real_t), intent(inout) :: Sarr(:,:)
        real(real_t), intent(out) :: Sout, Serr
        real(real_t) :: ak2, ak3, ak4, ak5, ak6, Stmp

        Stmp = S + b21 * h * dSdt
        call this%derivs(idx, t + a2 * h, Stmp, Sarr, ak2)
        Stmp = S + h * (b31*dSdt + b32*ak2)

        call this%derivs(idx, t + a3 * h, Stmp, Sarr, ak3)
        Stmp = S + h * (b41*dSdt + b42*ak2 + b43*ak3)

        call this%derivs(idx, t + a4 * h, Stmp, Sarr, ak4)
        Stmp = S + h * (b51*dSdt + b52*ak2 + b53*ak3 + b54*ak4)

        call this%derivs(idx, t + a5 * h, Stmp, Sarr, ak5)
        Stmp = S + h * (b61*dSdt + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)

        call this%derivs(idx, t + a6 * h, Stmp, Sarr, ak6)
        Sout = S + h * (c1*dSdt + c3*ak3 + c4*ak4 + c6*ak6)
        Serr = h * (dc1*dSdt + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)
    end subroutine rkck_sub

    subroutine rkqs_sub(this, idx, S, Sarr, dSdt, t, htry, hdid, hnext)
        class(runge_kutta_t) :: this
        integer, intent(in) :: idx
        real(real_t), intent(inout) :: t, S, Sarr(:,:)
        real(real_t), intent(in) :: dSdt, htry 
        real(real_t), intent(out) :: hdid, hnext
        real(real_t) :: errmax, h, htmp, tnew, Serr, Stmp
        real(real_t), parameter :: safety = 0.9_real_t, pgrow = -0.2_real_t 
        real(real_t), parameter :: pshrnk = -0.25_real_t
        real(real_t), parameter :: errcon = 1.89e-4_real_t
        real(real_t), parameter :: tiny_real_t = 1.e-30_real_t

        h = htry
        do
            call rkck_sub(this, idx, S, Sarr, dSdt, t, h, Stmp, Serr)
            errmax = abs(Serr/this%yscal) / this%eps
            if (errmax <= 1.0_real_t) exit

            htmp = safety * h * errmax**pshrnk
            h = sign( max(abs(htmp), 0.1_real_t*abs(h)), h )
            tnew = t + h

            if ( abs(tnew - t) < tiny_real_t ) then
                error stop 'stepsize underflow in rkqs'
            end if 
        end do

        if (errmax > errcon) then
            hnext = safety * h * errmax**pgrow
        else
            hnext = 5.0_real_t * h
        end if

        hdid = h; t = t + h; S = Stmp
    end subroutine rkqs_sub

    subroutine rkm_sub(this, idx, S, Sarr, dSdt, t, htry, hdid, hnext)
        class(runge_kutta_t) :: this
        integer, intent(inout) :: idx
        real(real_t), intent(inout) :: t, S, Sarr(:,:)
        real(real_t), intent(in) :: htry, dSdt
        real(real_t), intent(out) :: hdid, hnext

        if (idx < 1) then 
            error stop "Index has to be equal or greater than one."
        else if (idx > 9999) then 
            error stop "Index greater than the max value possible."
        end if 

        Sarr(1:4, 1:idx - 1) = Sarr_global(1:4, 1:idx - 1)
        call rkqs_sub(this, idx, S, Sarr, dSdt, t, htry, hdid, hnext)
        Sarr_global(1:4, idx) = [hnext, t, S, dSdt]
    end subroutine rkm_sub
   
end module runge_kutta_m
