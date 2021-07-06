subroutine heat_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,lambda, method,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn, method
    DOUBLE PRECISION xlower,ylower,dx,dy, lambda
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m, blockno
    DOUBLE PRECISION xc,yc

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
!!            xc = xlower + (i-0.5)*dx
!!            yc = ylower + (j-0.5)*dy
            rhs(i,j,1) = lambda*q(i,j,1)
        end do
    end do

end subroutine heat_fort_rhs


subroutine heat_update_q(mbc,mx,my,meqn,mfields,rhs,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j

    do i = 1,mx
        do j = 1,my
            q(i,j,1) = rhs(i,j,1)
        end do
    end do
100 format(2F16.8)

end subroutine heat_update_q

