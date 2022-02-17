subroutine poisson_fort_rhs(blockno, mbc,mx,my,mfields, & 
                           xlower,ylower,dx,dy,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision q

    INTEGER i,j, m
    DOUBLE PRECISION xc,yc, poisson_qexact_rhs
    INTEGER blockno

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            q = poisson_qexact_rhs(xc,yc)
            do m = 1,mfields
                rhs(i,j,m) =  q
            end do
        end do
    end do

end subroutine poisson_fort_rhs


