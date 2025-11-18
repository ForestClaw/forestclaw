DOUBLE PRECISION function poisson_qexact(x,y,z)
    IMPLICIT NONE

    DOUBLE PRECISION x,y,z
    
    INTEGER example
    COMMON /comm_example/ example

    INTEGER flag
    DOUBLE PRECISION grad(3), q, qlap


    flag = 0  !! Don't compute the gradient
    call poisson_qexact_complete(example,x,y,z,q,qlap,grad,flag)

    poisson_qexact = q

end function poisson_qexact

DOUBLE PRECISION function poisson_qexact_rhs(x,y,z)
    implicit none

    double precision x,y,z

    INTEGER example
    COMMON /comm_example/ example

    integer flag
    double precision q,qlap,b, grad_q(3), grad_beta(3)

    CALL poisson_fort_beta(x,y,z,b,grad_beta)

    flag = 2
    CALL poisson_qexact_complete(example,x,y,z,q,qlap,grad_q,flag)

    poisson_qexact_rhs = grad_beta(1)*grad_q(1) +  grad_beta(2)*grad_q(2) + grad_beta(3)*grad_q(3) + b*qlap

END FUNCTION poisson_qexact_rhs



SUBROUTINE poisson_qexact_gradient(x,y,z,q,grad)
    IMPLICIT NONE

    DOUBLE PRECISION x,y,z, q, grad(3)

    INTEGER flag
    DOUBLE PRECISION qlap

    INTEGER example
    COMMON /comm_example/ example

    flag = 1
    CALL poisson_qexact_complete(example,x,y,z,q,qlap,grad,flag)

END SUBROUTINE poisson_qexact_gradient



SUBROUTINE poisson_qexact_complete(example,x,y,z,q,qlap,grad,flag)
    use hsmooth_mod, only : m_polar, x0_polar, y0_polar
    IMPLICIT NONE

    DOUBLE PRECISION x,y,z, q, qlap, grad(3)
    INTEGER flag, example

!!    INTEGER example
!!    COMMON /comm_example/ example

    DOUBLE PRECISION alpha,x0,y0,z0,a,b,c
    COMMON /comm_rhs/ alpha,x0,y0,z0,a,b,c

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION r, r2, theta
    double precision hsmooth, h_grad(3), hsmooth_laplacian
    DOUBLE PRECISION qx,qy,qz, dqdr, t1(3), t2(3)
    double precision q1, qx1, qy1, qlap1, x0p,y0p
    integer id

    if (example .eq. 0) then
        q = x**2 + y**2 + z**2
        qx = 2*x
        qy = 2*y
        qz = 2*z
        qlap = 6
    elseif (example .eq. 1) then
        !! example in polar coordinates (r)
        r2 = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
        q1 = exp(-alpha/2.d0*r2)
        q = q1 + 1
        if (flag .ge. 1) then
            r = sqrt(r2)
            t1(1) = (x-x0)/r
            t1(2) = (y-y0)/r
            t1(3) = (z-z0)/r
            dqdr = -alpha*r*q1
            qx = dqdr*t1(1)  !! Cartesian components of gradient
            qy = dqdr*t1(2)
            qz = dqdr*t1(3)
            if (flag .eq. 2) then
                qlap = alpha*exp(-alpha/2.d0*r2)*(alpha*r2 - 2)
            endif
        endif
    elseif (example .eq. 2) then
        !! Example in Cartesian coordinates
        q = cos(pi*a*x)*cos(pi*b*y)*cos(pi*c*z)
        if (flag .ge. 1) then
            qx = -pi*a*sin(pi*a*x)*cos(pi*b*y)*cos(pi*c*z)
            qy = -pi*b*cos(pi*a*x)*sin(pi*b*y)*cos(pi*c*z)
            qz = -pi*c*sin(pi*a*x)*cos(pi*b*y)*sin(pi*c*z)
            if (flag .eq. 2 ) then
                qlap = -pi**2*(a**2 + b**2 + c**2)*cos(pi*a*x)*cos(pi*b*y)*cos(pi*c*z)
            endif
        endif
    elseif (example .eq. 3) then
        !! d/dx = (y - 1) y (z - 1) z (-e^(x y z)) (x^2 y z + x (2 - y z) - 1)
        !! d/dy = (x - 1) x (z - 1) z (-e^(x y z)) (x y^2 z + y (2 - x z) - 1)
        !! d/dz = (x - 1) x (y - 1) y (-e^(x y z)) (x y z^2 + z (2 - x y) - 1)
        q = (1-x)*x*(1-y)*y*(1-z)*z*exp(x*y*z)

        if (flag .ge. 1) then
            qx = -(y-1)*y*(z-1)*z*exp(x*y*z)*(x**2*y*z + x*(2-y*z) - 1)
            qy = -(x-1)*x*(z-1)*z*exp(x*y*z)*(x*y**2*z + y*(2-x*z) - 1)
            qz = -(x-1)*x*(y-1)*y*exp(x*y*z)*(x*y*z**2 + z*(2-x*y) - 1)
            if (flag .eq. 2) then
                !! Δ((1 - x) x (1 - y) y (1 - z) z exp(x y z)) 
                !! = e^(x y z) * (
                !!     -(x - 1) x^3 (y - 1) y^3 (z - 1) z 
                !!     - (x - 1) x^3 (y - 1) y (z - 1) z^3 
                !!     - 2 (x - 1) x^2 (y - 1) y^2 (z - 1) 
                !!     - 2 (x - 1) x^2 (y - 1) y^2 z 
                !!     - 2 (x - 1) x^2 y (z - 1) z^2 
                !!     - 2 (x - 1) x^2 (y - 1) (z - 1) z^2 
                !!     - ((x - 1) x (y - 1) y^3 (z - 1) z^3) 
                !!     - 2 (x - 1) (y - 1) y^2 (z - 1) z^2 
                !!     - 2 x (y - 1) y^2 (z - 1) z^2 
                !!     - 2 (x - 1) x (y - 1) y 
                !!     - 2 (x - 1) x (z - 1) z 
                !!     - 2 (y - 1) y (z - 1) z
                !! )
                qlap = exp(x*y*z)*( &
                        -(x-1)*x**3*(y-1)*y**3*(z-1)*z &
                        -(x-1)*x**3*(y-1)*y*(z-1)*z**3 &
                        -2*(x-1)*x**2*(y-1)*y**2*(z-1) & 
                        -2*(x-1)*x**2*(y-1)*y**2*z &
                        -2*(x-1)*x**2*y*(z-1)*z**2 &
                        -2*(x-1)*x**2*(y-1)*(z-1)*z**2 &
                        -((x-1)*x*(y-1)*y**3*(z-1)*z**3) &
                        -2*(x-1)*(y-1)*y**2*(z-1)*z**2 &
                        -2*x*(y-1)*y**2*(z-1)*z**2 &
                        -2*(x-1)*x*(y-1)*y &
                        -2*(x-1)*x*(z-1)*z &
                        -2*(y-1)*y*(z-1)*z &
                      )
            endif
        endif
    elseif (example .eq. 4) then
        q = 0
        qx = 0
        qy = 0
        qz = 0
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
        grad(3) = qz
    endif

end subroutine poisson_qexact_complete

double precision function sech(x)
    implicit none

    double precision x

    sech = 1.d0/cosh(x)
end function sech

subroutine poisson_fort_beta(x,y,z,b,grad)
    implicit none

    double precision x,y,z,b,grad(3)

    integer beta_choice
    common /comm_beta/ beta_choice

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION bx, by, bz

    if (beta_choice .eq. 0) then
        b = 1
        bx = 0
        by = 0
        bz = 0
    elseif (beta_choice .eq. 1) then
        b = cos(pi*x)*cos(pi*y)*cos(pi*z) + 2
        bx = -pi*sin(pi*x)*cos(pi*y)*cos(pi*z)
        by = -pi*cos(pi*x)*sin(pi*y)*cos(pi*z)
        bz = -pi*cos(pi*x)*cos(pi*y)*sin(pi*z)
    elseif (beta_choice .eq. 2) then
        b = 1 + x*y*z
        bx = y*z
        by = x*z
        bz = x*y
    endif

    grad(1) = bx
    grad(2) = by
    grad(3) = bz

end subroutine poisson_fort_beta

