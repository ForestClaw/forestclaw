subroutine allencahn_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,dt, method,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn, method
    DOUBLE PRECISION xlower,ylower,dx,dy, dt
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m, blockno
    DOUBLE PRECISION xc,yc, lambda, D, u

    lambda = -1.d0/dt

    D = 0.001d0

!!do j = 1,my
!!         do i = 1,mx
!!            u = q(i,j,1)
!!
!!c           # Treat only diffusion term explicitly
!!            rhs(i,j,1) = (u-u**3)/D**2
!!         enddo
!!      enddo

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            u = q(i,j,1)
            rhs(i,j,1) = lambda*u - (u-u**3)/D**2
        end do
    end do

end subroutine allencahn_fort_rhs


subroutine allencahn_update_q(mbc,mx,my,meqn,mfields,rhs,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            q(i,j,1) = rhs(i,j,1)
        end do
    end do
100 format(2F16.8)

end subroutine allencahn_update_q

