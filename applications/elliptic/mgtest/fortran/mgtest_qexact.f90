double precision function mgtest_qexact(x,y)
    IMPLICIT NONE

    DOUBLE PRECISION x,y

    INTEGER rhs_choice
    COMMON /common_rhs_int/ rhs_choice

    DOUBLE PRECISION alpha,x0,y0,a,b
    COMMON /common_rhs_double/ alpha,x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /common_pi/ pi, pi2

    INTEGER i,j
    DOUBLE PRECISION q, r, r2, r0, hsmooth, hsmooth_deriv


    if (rhs_choice .eq. 1) then
        r2 = (x-x0)**2 + (y-y0)**2
        q = exp(-alpha/2.d0*r2)
    elseif (rhs_choice .eq. 2) then
        q = sin(pi*a*x)*sin(pi*b*y)
    elseif (rhs_choice .eq. 3) then
        write(6,*) 'qexact.f : Exact solution not yet implemented for  choice == 3'
        stop
        r0 = 0.25
        r = sqrt((x-0.5)**2 + (y-0.5)**2)
!!      q(i,j) = Hsmooth(r + r0) - Hsmooth(r - r0)
        q = hsmooth_deriv(r + r0) - hsmooth_deriv(r - r0)
    endif

    mgtest_qexact = q

end function mgtest_qexact