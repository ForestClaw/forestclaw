subroutine mgtest_setprob(rhs_in, x0_in, y0_in, a_in, b_in)
    implicit none


    INTEGER rhs_in
    DOUBLE PRECISION x0_in, y0_in, a_in, b_in

    INTEGER rhs_choice
    COMMON /common_rhs_int/ rhs_choice

    DOUBLE PRECISION x0,y0,a,b
    COMMON /common_rhs_double/ x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /common_pi/ pi, pi2

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    rhs_choice = rhs_in

    x0 = x0_in
    y0 = y0_in

    a = a_in
    b = b_in

end subroutine mgtest_setprob