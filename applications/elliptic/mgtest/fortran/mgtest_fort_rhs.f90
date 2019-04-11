subroutine mgtest_fort_rhs(mbc,mx,my,xlower,ylower,dx,dy,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER rhs_choice
    COMMON /common_rhs_int/ rhs_choice

    DOUBLE PRECISION x0,y0,a,b
    COMMON /common_rhs_double/ x0,y0,a,b

    DOUBLE PRECISION pi,pi2
    COMMON /common_pi/ pi, pi2

    INTEGER i,j
    DOUBLE PRECISION x,y,r2

    do i = 1,mx
        do j = 1,my
            x = xlower + (i-0.5)*dx
            y = ylower + (j-0.5)*dy
            if (rhs_choice .eq. 1) then
                r2 = (x-x0)**2 + (y-y0)**2
                q(i,j) = exp(-20*r2)
            else
                q(i,j) = cos(pi*a*x)*cos(pi*b*y)
            endif
        end do
    end do

end subroutine mgtest_fort_rhs



