SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    IMPLICIT NONE

    INTEGER meqn,mbc,mx,my,maux
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

    INTEGER i,j,mq
    DOUBLE PRECISION xlow,ylow, wl
    integer blockno, fc2d_clawpack46_get_block

    blockno = fc2d_clawpack46_get_block()

    do mq = 1,meqn
        DO  i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            DO j = 1-mbc,my+mbc
                ylow = ylower + (j-1)*dy
                call cellave2(blockno,xlow,ylow,dx,dy,wl)
                q(1,i,j) = wl
            END DO
        END DO
    END DO

END SUBROUTINE clawpack5_qinit
