subroutine mgtest_fort_rhs(blockno, mbc,mx,my,xlower,ylower,dx,dy,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER rhs_choice
    COMMON /common_rhs_int/ rhs_choice

    DOUBLE PRECISION alpha,x0,y0,a,b
    COMMON /common_rhs_double/ alpha,x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /common_pi/ pi, pi2

    INTEGER*8 cont, get_context


    INTEGER i,j
    DOUBLE PRECISION xc,yc, xc1, yc1, zc1, r, r2, r0
    DOUBLE PRECISION hsmooth, hsmooth_deriv
    integer fclaw2d_blockno, blockno

    cont = get_context()

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            if (rhs_choice .eq. 1) then
                r2 = (xc1-x0)**2 + (yc1-y0)**2                
                q(i,j) = exp(-alpha/2.d0*r2)*(alpha**2*r2 - 2*alpha)
            elseif (rhs_choice .eq. 2) then
!!               # q(x,y) = sin(pi*a*x)*sin(pi*b*y) 
                q(i,j) = -(pi**2*(a**2 + b**2))*sin(pi*a*xc1)*sin(pi*b*yc1)               
            elseif (rhs_choice .eq. 3) then
                r0 = 0.25
                r = sqrt((xc1-0.5)**2 + (yc1-0.5)**2)
!!                q(i,j) = Hsmooth(r + r0) - Hsmooth(r - r0)
                q(i,j) = hsmooth_deriv(r + r0) - hsmooth_deriv(r - r0)
            endif
        end do
    end do

end subroutine mgtest_fort_rhs


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


