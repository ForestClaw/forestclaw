DOUBLE PRECISION function mgtest_qexact(x,y)
    IMPLICIT NONE

    DOUBLE PRECISION x,y
    
    INTEGER flag
    DOUBLE PRECISION grad(2), q, qlap


    flag = 0  !! Don't compute the gradient
    call mgtest_qexact_complete(x,y,q,qlap,grad,flag)

    mgtest_qexact = q

end function mgtest_qexact

DOUBLE PRECISION function mgtest_qexact_rhs(x,y)
    implicit none

    double precision x,y

    integer flag
    double precision q,qlap,grad(2)

    flag = 2
    CALL mgtest_qexact_complete(x,y,q,qlap,grad,flag)

    mgtest_qexact_rhs = qlap

END FUNCTION mgtest_qexact_rhs


SUBROUTINE mgtest_qexact_gradient(x,y,q,grad)
    IMPLICIT NONE

    DOUBLE PRECISION x,y, q, grad(2)

    INTEGER flag
    DOUBLE PRECISION qlap

    flag = 1
    CALL mgtest_qexact_complete(x,y,q,qlap,grad,flag)

END SUBROUTINE mgtest_qexact_gradient



SUBROUTINE mgtest_qexact_complete(x,y,q,qlap,grad,flag)
    IMPLICIT NONE

    DOUBLE PRECISION x,y, q, qlap, grad(2)
    INTEGER flag

    INTEGER example
    COMMON /comm_example/ example

    DOUBLE PRECISION alpha,x0,y0,a,b
    COMMON /comm_rhs/ alpha,x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    INTEGER i,j
    DOUBLE PRECISION r, r2, r0, hsmooth, hsmooth_deriv
    DOUBLE PRECISION rx,ry,qx,qy

    if (example .eq. 1) then
        r2 = (x-x0)**2 + (y-y0)**2
        q = exp(-alpha/2.d0*r2)
        if (flag .eq. 1) then
            rx = (x-x0)/sqrt(r2)
            ry = (y-y0)/sqrt(r2)
            qx = -alpha*r*rx*q
            qy = -alpha*r*ry*q
        endif
        if (flag .eq. 2) then
             qlap = exp(-alpha/2.d0*r2)*(alpha**2*r2 - 2*alpha)
        endif
    elseif (example .eq. 2) then
        q = sin(pi*a*x)*sin(pi*b*y)
        if (flag .eq. 1) then
            qx = pi*a*cos(pi*a*x)*sin(pi*b*y)
            qy = pi*b*sin(pi*a*x)*cos(pi*b*y)
        endif
        if (flag .eq. 2 ) then
            qlap = -(pi**2*(a**2 + b**2))*sin(pi*a*x)*sin(pi*b*y)
        endif
    elseif (example .eq. 3) then
        !! Variable coefficient problem
    elseif (example .eq. 4) then
        write(6,*) 'mgtest_qexact.f : Exact solution not yet implemented for  choice == 3'
        stop
        r0 = 0.25
        r = sqrt((x-0.5)**2 + (y-0.5)**2)
!!      q(i,j) = Hsmooth(r + r0) - Hsmooth(r - r0)
        q = hsmooth_deriv(r + r0) - hsmooth_deriv(r - r0)
    endif

    if (flag .eq. 1) then
        grad(1) = qx
        grad(2) = qy
    endif

end subroutine mgtest_qexact_complete


double precision function Hsmooth(r)
    implicit none

    double precision r, a

    a = 0.015625d0
    Hsmooth = (tanh(r/a) + 1)/2.

end function Hsmooth

double precision function Hsmooth_deriv(r)
    implicit none

    double precision r, a

    a = 0.015625d0
    Hsmooth_deriv = (1/a)*(1./cosh(r/a))**2/2.

end function Hsmooth_deriv


