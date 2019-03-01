module mymod
    use mpmod
    implicit none

    real(mp), parameter :: pi = 2.0_mp*asin(1.0_mp)
    real(mp), parameter :: par = 0.02654_mp
    real(mp), parameter :: k_b = 1.38e-16_mp   !постоянная Больцмана
    real(mp), parameter :: aem = 1.66e-24_mp   !атомная единица массы
    real(mp), parameter :: h   = 6.63e-27_mp
    real(mp), parameter :: m_e = 9.11e-28_mp
    real(mp), parameter :: q_e = -4.80e-10_mp
    real(mp), parameter :: v_c =    2.998e10_mp
    real(mp) :: imp_r

    real(mp) :: v_esc, r_0, R_s, r_end
end module mymod
