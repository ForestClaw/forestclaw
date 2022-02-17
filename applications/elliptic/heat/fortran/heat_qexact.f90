DOUBLE PRECISION function heat_qexact(x,y)
    IMPLICIT NONE

    DOUBLE PRECISION x,y
    
    INTEGER flag
    DOUBLE PRECISION grad(2), q, qlap


    flag = 0  !! Don't compute the gradient
    call heat_qexact_complete(x,y,q,qlap,grad,flag)

    heat_qexact = q

end function heat_qexact

DOUBLE PRECISION function heat_qexact_rhs(x,y)
    implicit none

    double precision x,y

    integer flag
    double precision q,qlap,b, grad_q(2), grad_beta(2)

    CALL heat_fort_beta(x,y,b,grad_beta)

    flag = 2
    CALL heat_qexact_complete(x,y,q,qlap,grad_q,flag)

    heat_qexact_rhs = grad_beta(1)*grad_q(1) +  grad_beta(2)*grad_q(2) + b*qlap

END FUNCTION heat_qexact_rhs



SUBROUTINE heat_qexact_gradient(x,y,q,grad)
    IMPLICIT NONE

    DOUBLE PRECISION x,y, q, grad(2)

    INTEGER flag
    DOUBLE PRECISION qlap

    flag = 1
    CALL heat_qexact_complete(x,y,q,qlap,grad,flag)

END SUBROUTINE heat_qexact_gradient



SUBROUTINE heat_qexact_complete(x,y,q,qlap,grad,flag)
    use hsmooth_mod, only : m_polar, x0_polar, y0_polar
    IMPLICIT NONE

    DOUBLE PRECISION x,y, q, qlap, grad(2)
    INTEGER flag

    INTEGER example
    COMMON /comm_example/ example

    DOUBLE PRECISION alpha,x0,y0,a,b
    COMMON /comm_rhs/ alpha,x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION r, r2, theta
    double precision hsmooth, h_grad(2), hsmooth_laplacian
    DOUBLE PRECISION qx,qy, dqdr, t1(2), t2(2)
    double precision q1, qx1, qy1, qlap1, x0p,y0p
    integer id

    if (example .eq. 0) then
        q = x
        qx = 1
        qy = 1
        qlap = 0
    elseif (example .eq. 1) then
        !! example in polar coordinates (r)
        r2 = (x-x0)**2 + (y-y0)**2
        q1 = exp(-alpha/2.d0*r2)
        q = q1 + 1
        if (flag .ge. 1) then
            r = sqrt(r2)
            t1(1) = (x-x0)/r
            t1(2) = (y-y0)/r
            dqdr = -alpha*r*q1
            qx = dqdr*t1(1)  !! Cartesian components of gradient
            qy = dqdr*t1(2)
            if (flag .eq. 2) then
                qlap = alpha*exp(-alpha/2.d0*r2)*(alpha*r2 - 2)
            endif
        endif
    elseif (example .eq. 2) then
        !! Example in Cartesian coordinates
        q = cos(pi*a*x)*cos(pi*b*y)
        if (flag .ge. 1) then
            qx =  pi*a*cos(pi*a*x)*cos(pi*b*y)
            qy = -pi*b*sin(pi*a*x)*sin(pi*b*y)
            if (flag .eq. 2 ) then
                qlap = -(pi**2*(a**2 + b**2))*cos(pi*a*x)*cos(pi*b*y)
            endif
        endif
    elseif (example .eq. 3) then
        !! d/(dx)((1 - x) x (1 - y) y exp(x y)) = (y - 1) y e^(x y) (x^2 y - (y - 2) x - 1)
        !! d/(dy)((1 - x) x (1 - y) y exp(x y)) = (x - 1) x e^(x y) (x y^2 - (x - 2) y - 1)
        q = (1-x)*x*(1-y)*y*exp(x*y)

        if (flag .ge. 1) then
            qx = (y - 1)*y*exp(x*y)*(x**2*y - (y - 2)*x - 1)
            qy = (x - 1)*x*exp(x*y)*(x*y**2 - (x - 2)*y - 1)
            if (flag .eq. 2) then
                !!Δ((1 - x) x (1 - y) y e^(x y)) = 
                !! e^(x y) (x^4 (y - 1) y 
                !! - x^3 (y^2 - 5 y + 2) 
                !! + x^2 (y^4 - y^3 - 4 y + 4) 
                !! - x (y^4 - 5 y^3 + 4 y^2 + 2) 
                !! - 2 (y - 1)^2 y)
                qlap = exp(x*y)*(x**4*(y - 1)*y & 
                                 - x**3*(y**2 - 5*y + 2) & 
                                 + x**2*(y**4 - y**3 - 4*y + 4) & 
                                 - x*(y**4 - 5*y**3 + 4*y**2 + 2) &
                                 - 2*(y - 1)**2*y)
            endif
        endif
    elseif (example .eq. 4) then
        q = 0
        qx = 0
        qy = 0
        qlap = 0
        do id = 1,m_polar
            x0p = x0_polar(id)
            y0p = y0_polar(id)
            r = sqrt((x-x0p)**2 + (y-y0p)**2)
            theta = atan2(y-y0p,x-x0p)        
            q1 = 1 - hsmooth(id,r,theta)
            if (flag .ge. 1) then
                !! Assume mapping is T(r,theta)
                t1(1) = cos(theta)
                t1(2) = sin(theta) 
                t2(1) = -r*sin(theta)
                t2(2) = r*cos(theta)

                !! Cartesian components of the gradient
                call hsmooth_grad(id,r,theta,h_grad)
                qx1 = -(h_grad(1)*t1(1) + h_grad(2)*t2(1))   
                qy1 = -(h_grad(1)*t1(2) + h_grad(2)*t2(2))   

                !! Laplacian
                if (flag .eq. 2) then
                    qlap1 = -hsmooth_laplacian(id,r,theta)
                endif
            endif
            q = q + q1
            qx = qx + qx1
            qy = qy + qy1
            qlap = qlap + qlap1
        enddo
    endif
    if (flag .ge. 1) then
        grad(1) = qx
        grad(2) = qy
    endif

end subroutine heat_qexact_complete

double precision function sech(x)
    implicit none

    double precision x

    sech = 1.d0/cosh(x)
end function sech

subroutine heat_fort_beta(x,y,b,grad)
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

end subroutine heat_fort_beta

