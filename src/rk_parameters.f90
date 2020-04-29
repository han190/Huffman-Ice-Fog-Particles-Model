module rk_parameters_m
    use real_m
    implicit none

    real(real_t), parameter :: a2  = 0.2_real_t
    real(real_t), parameter :: a3  = 0.3_real_t
    real(real_t), parameter :: a4  = 0.6_real_t
    real(real_t), parameter :: a5  = 1.0_real_t
    real(real_t), parameter :: a6  = 0.875_real_t
    real(real_t), parameter :: b21 = 0.2_real_t
    real(real_t), parameter :: b31 = 3.0_real_t / 40.0_real_t
    real(real_t), parameter :: b32 = 9.0_real_t / 40.0_real_t
    real(real_t), parameter :: b41 = 0.3_real_t
    real(real_t), parameter :: b42 = -0.9_real_t
    real(real_t), parameter :: b43 = 1.2_real_t
    real(real_t), parameter :: b51 = -11.0_real_t / 54.0_real_t
    real(real_t), parameter :: b52 = 2.5_real_t
    real(real_t), parameter :: b53 = -70.0_real_t / 27.0_real_t
    real(real_t), parameter :: b54 = 35.0_real_t / 27.0_real_t
    real(real_t), parameter :: b61 = 1631.0_real_t / 55296.0_real_t
    real(real_t), parameter :: b62 = 175.0_real_t / 512.0_real_t
    real(real_t), parameter :: b63 = 575.0_real_t / 13824.0_real_t
    real(real_t), parameter :: b64 = 44275.0_real_t / 110592.0_real_t
    real(real_t), parameter :: b65 = 253.0_real_t / 4096.0_real_t
    real(real_t), parameter :: c1  = 37.0_real_t / 378.0_real_t
    real(real_t), parameter :: c3  = 250.0_real_t / 621.0_real_t
    real(real_t), parameter :: c4  = 125.0_real_t / 594.0_real_t
    real(real_t), parameter :: c6  = 512.0_real_t / 1771.0_real_t
    real(real_t), parameter :: dc1 = c1 - 2825.0_real_t / 27648.0_real_t
    real(real_t), parameter :: dc3 = c3 - 18575.0_real_t / 48384.0_real_t
    real(real_t), parameter :: dc4 = c4 - 13525.0_real_t / 55296.0_real_t
    real(real_t), parameter :: dc5 = -277.0_real_t / 14336.0_real_t
    real(real_t), parameter :: dc6 = c6 - 0.25_real_t

end module rk_parameters_m