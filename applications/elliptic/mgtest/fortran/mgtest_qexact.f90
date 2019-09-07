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
    double precision q,qlap,b, grad_q(2), grad_beta(2)

    CALL mgtest_beta(x,y,b,grad_beta)

    flag = 2
    CALL mgtest_qexact_complete(x,y,q,qlap,grad_q,flag)

    mgtest_qexact_rhs = grad_beta(1)*grad_q(1) +  grad_beta(2)*grad_q(2) + b*qlap

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
    DOUBLE PRECISION r, r2, r0
    double precision hsmooth, hsmooth_deriv, hsmooth_laplacian
    DOUBLE PRECISION rx,ry,qx,qy, grad_beta

    if (example .eq. 0) then
        q = x + y
        qx = 1
        qy = 1
        qlap = 0
    elseif (example .eq. 1) then
        r2 = (x-x0)**2 + (y-y0)**2
        q = exp(-alpha/2.d0*r2) + 1
        if (flag .ge. 1) then
            rx = (x-x0)/sqrt(r2)
            ry = (y-y0)/sqrt(r2)
            qx = -alpha*r*rx*q
            qy = -alpha*r*ry*q
            if (flag .eq. 2) then
                qlap = exp(-alpha/2.d0*r2)*(alpha**2*r2 - 2*alpha)
            endif
        endif
    elseif (example .eq. 2) then
        q = sin(pi*a*x)*cos(pi*b*y)
        if (flag .ge. 1) then
            qx = pi*a*cos(pi*a*x)*cos(pi*b*y)
            qy = -pi*b*sin(pi*a*x)*sin(pi*b*y)
            if (flag .eq. 2 ) then
                qlap = -(pi**2*(a**2 + b**2))*sin(pi*a*x)*cos(pi*b*y)
            endif
        endif
    elseif (example .eq. 3) then
        !! Variable coefficient problem
        !! d/(dy)((1 - x) x (1 - y) y exp(x y)) = (x - 1) x e^(x y) (x y^2 - (x - 2) y - 1)
        q = (1-x)*x*(1-y)*y*exp(x*y)

        if (flag .ge. 1) then
            qx = (x - 1)*x*exp(x*y)*(x*y**2 - (x - 2)*y - 1)
            qy = (y - 1)*x*exp(x*y)*(x*y**2 - (x - 2)*y - 1)
            if (flag .eq. 2) then
                !!Δ((1 - x) x (1 - y) y e^(x y)) = e^(x y) (x^4 (y - 1) y - 
                !!x^3 (y^2 - 5 y + 2) +  x^2 
                !!(y^4 - y^3 - 4 y + 4) - x (y^4 - 5 y^3 + 4 y^2 + 2) - 2 (y - 1)^2 y)
                qlap = exp(x*y)*(x**4*(y - 1)*y & 
                                 - x**3*(y**2 - 5*y + 2) & 
                                 + x**2*(y**4 - y**3 - 4*y + 4) & 
                                 - x*(y**4 - 5*y**3 + 4*y**2 + 2) &
                                 - 2*(y - 1)**2*y)
            endif
        endif
    elseif (example .eq. 4) then
        r0 = 0.25
        r = sqrt((x-0.5)**2 + (y-0.5)**2)
        q = hsmooth(r + x0) - hsmooth(r - y0)
        if (flag .ge. 1) then
            rx = (x-x0)/sqrt(r2)
            ry = (y-y0)/sqrt(r2)
            qx = (hsmooth_deriv(r + r0) - hsmooth_deriv(r - r0))*rx
            qy = (hsmooth_deriv(r + r0) - hsmooth_deriv(r - r0))*ry
            if (flag .eq. 2) then
                qlap = hsmooth_laplacian(r+r0) - hsmooth_laplacian(r-r0)
            endif
        endif
    endif
    if (flag .ge. 1) then
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

!! Compute derivatives with respect to r
DOUBLE PRECISION function Hsmooth_deriv(r)
    implicit none

    double precision r,a, sech2

    a = 0.015625d0
    sech2 = (1./cosh(r/a))**2
    Hsmooth_deriv = sech2/(2*a)

end function Hsmooth_deriv

!! Compute Laplacian 
double precision function Hsmooth_laplacian(r)
    implicit none

    double precision r, a, sech2, a2

    a = 0.015625d0
    sech2 = (1.d0/cosh(r/a))**2
    a2 = a**2

    if (r .eq. 0) then
        write(6,*) 'r == 0)'
        stop
    endif
    Hsmooth_laplacian = (a-2*r*tanh(r/a))*sech2/(2*a2*r)

end function Hsmooth_laplacian


subroutine mgtest_beta(x,y,b,grad)
    implicit none

    double precision x,y,b,grad(2)

    integer beta_choice
    common /comm_beta/ beta_choice

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2
    
    DOUBLE PRECISION bx, by

    if (beta_choice .eq. 0) then
        b = 1
        bx = 0
        by = 0
    elseif (beta_choice .eq. 1) then
        b = cos(pi*x)*cos(pi*y) + 2
        bx = -pi*sin(pi*x)*cos(pi*y)
        by = -pi*cos(pi*x)*sin(pi*y)
    elseif (beta_choice .eq. 2) then
        b = 1 + x*y
        bx = y
        by = x
    endif

    grad(1) = bx
    grad(2) = by

end subroutine mgtest_beta

