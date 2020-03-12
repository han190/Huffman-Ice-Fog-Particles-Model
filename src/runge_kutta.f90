module runge_kutta_m 
    use real_m 
    implicit none 
    private

    integer, parameter :: nmax = 10000

    type, abstract, public :: runge_kutta_t 
        real(real_t) :: eps 
        real(real_t), pointer :: yscal(:)
    contains 
        procedure :: rk_init => initialize_sub 
        procedure :: rkck => rkck_sub 
        procedure :: rkqs => rkqs_sub 
        procedure(derivsx), deferred :: derivs 
    end type

    abstract interface

        subroutine derivsx(this, x, y, yarr, dydx)
            use real_m
            import :: runge_kutta_t
            implicit none

            class(runge_kutta_t) :: this
            real(real_t), intent(in) :: x
            real(real_t), intent(in) :: y(:)
            real(real_t), intent(inout) :: yarr(:,:)
            real(real_t), intent(out) :: dydx(:)
        end subroutine derivsx

    end interface

contains 

    subroutine initialize_sub(this, eps, yscal)
        class(runge_kutta_t) :: this 
        real(real_t), intent(in) :: eps 
        real(real_t), target, intent(in) :: yscal(:)

        this%eps = eps 
        this%yscal => yscal 
    end subroutine initialize_sub 

    subroutine rkck_sub(this, y, yarr, dydx, x, h, yout, yerr)
        use rk_parameters_m
        implicit none 

        class(runge_kutta_t) :: this
        real(real_t), intent(in) :: y(:), dydx(:), x, h
        real(real_t), intent(inout) :: yarr(:,:)
        real(real_t), intent(out) :: yout(:), yerr(:)
        real(real_t), dimension(size(y)) :: ak2, ak3, ak4, ak5, ak6, ytemp
        integer :: ndum
               
        ndum = size(y, dim = 1)
        ytemp = y + b21 * h * dydx
        call this%derivs(x + a2 * h, ytemp, yarr, ak2)
        ytemp = y + h * ( b31 * dydx + b32 * ak2 )

        call this%derivs(x + a3 * h, ytemp, yarr, ak3)
        ytemp = y + h * ( b41 * dydx + b42 * ak2 + b43 * ak3 )

        call this%derivs(x + a4 * h, ytemp, yarr, ak4)
        ytemp = y + h * ( b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4 )

        call this%derivs(x + a5 * h, ytemp, yarr, ak5)
        ytemp = y + h * ( &
            b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5 &
        )

        call this%derivs(x + a6 * h, ytemp, yarr, ak6)
        yout = y + h * ( &
            c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6 &
        )

        yerr = h * ( &
            dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6 &
        )
    end subroutine rkck_sub

    subroutine rkqs_sub(this, y, yarr, dydx, x, htry, hdid, hnext)
        class(runge_kutta_t) :: this
        real(real_t), intent(inout) :: y(:), yarr(:,:)
        real(real_t), intent(in) :: dydx(:)
        real(real_t), intent(inout) :: x
        real(real_t), intent(in) :: htry
        real(real_t), intent(out) :: hdid, hnext
        integer :: ndum
        real(real_t) :: errmax, h, htemp, xnew
        real(real_t), dimension(size(y)) :: yerr, ytemp
        real(real_t), parameter :: safety = 0.9_real_t, pgrow = -0.2_real_t 
        real(real_t), parameter :: pshrnk = -0.25_real_t
        real(real_t), parameter :: errcon = 1.89e-4_real_t
        real(real_t), parameter :: tiny_real_t = 1.e-30_real_t
        
        if ( .not. associated(this%yscal) ) then
            error stop ' yscal has not been associated.'
        end if 

        if ( size(this%yscal, dim=1) == size(y, dim=1) ) then 
            ndum = size(y, dim=1)
        else 
            error stop 'yscal and y have different dimensions'
        end if 

        h = htry
        do
            call rkck_sub(this, y, dydx, x, h, ytemp, yerr)
            errmax = maxval( abs( yerr(:) / this%yscal(:) ) ) / this%eps
            if (errmax <= 1.0_real_t) exit

            htemp = safety * h * errmax**pshrnk
            h = sign( max( abs(htemp), 0.1_real_t * abs(h) ), h )
            xnew = x + h

            if ( abs(xnew - x) < tiny_real_t ) then
                error stop 'stepsize underflow in rkqs'
            end if 
        end do

        if (errmax > errcon) then
            hnext = safety * h * errmax**pgrow
        else
            hnext = 5.0_real_t * h
        end if

        hdid = h; x = x + h; y(:) = ytemp(:)
    end subroutine rkqs_sub
   
end module runge_kutta_m