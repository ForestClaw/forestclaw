SUBROUTINE clawpack46_qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower, & 
        dx,dy,dz, q,maux,aux)
    IMPLICIT NONE

    INTEGER :: meqn, mbc, mx, my, mz, maux
    DOUBLE PRECISION :: xlower, ylower, zlower, dx, dy, dz
    DOUBLE PRECISION ::   q(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,maux)

    INTEGER :: i, j, k, mq
    DOUBLE PRECISION :: xlow, ylow, zlow, w   

    integer :: blockno, fc3d_clawpack46_get_block

    blockno = fc3d_clawpack46_get_block()

    DO mq = 1,meqn
        do k = 1-mbc,mz+mbc
            zlow = zlower + (k-1)*dz
            DO j = 1-mbc,my+mbc
                ylow = ylower + (j-1)*dy
                DO i = 1-mbc,mx+mbc
                    xlow = xlower + (i-1)*dx
                    CALL cellave3(blockno,xlow,ylow,zlow,dx,dy,dz,w)
                    q(i,j,k,1) = w
                end do                
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE clawpack46_qinit
