subroutine heat_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,dt, method,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn, method
    DOUBLE PRECISION xlower,ylower,dx,dy, dt
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m, blockno
    DOUBLE PRECISION xc,yc

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            if (method .eq. 1) then
                rhs(i,j,1) = q(i,j,1)
            else
                write(6,*) 'heat_fort_rhs : no valid RHS side specified'
                stop
            endif
        end do
    end do

end subroutine heat_fort_rhs


subroutine heat_update_q(mbc,mx,my,meqn,mfields,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m, blockno
    DOUBLE PRECISION tmp

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            q(i,j,1) = rhs(i,j,1)
        end do
    end do

end subroutine heat_update_q

