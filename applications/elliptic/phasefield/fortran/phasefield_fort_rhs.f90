subroutine phasefield_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,dt, method,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn, method
    DOUBLE PRECISION xlower,ylower,dx,dy, dt
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    DOUBLE PRECISION S, alpha, m_parm, xi, k, gamma
    COMMON /comm_parms/ S, alpha, m_parm, xi, k, gamma

    INTEGER i,j, m, blockno
    DOUBLE PRECISION lambda, u, phi, g0, g, s1, s3, beta

    lambda = -1.d0/dt
    beta = xi**2/m_parm;

    do j = 1,my
        do i = 1,mx
            u = q(i,j,1)
            phi = q(i,j,2)

            g0 = phi*(1-phi)
            g = g0**2

            s1 = 30*g/S
            s3 = g0*(phi-0.5d0)

            rhs(i,j,1) = lambda*(u + s1*phi)
            rhs(i,j,2) = lambda*beta*phi - s3
        end do
    end do

end subroutine phasefield_fort_rhs


subroutine phasefield_update_q(mbc,mx,my,meqn,mfields,rhs,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m

    do m = 1,2
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                q(i,j,1) = rhs(i,j,m)
            end do
        end do
    end do
100 format(2F16.8)

end subroutine phasefield_update_q

