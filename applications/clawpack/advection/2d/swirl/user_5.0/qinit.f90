SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    IMPLICIT NONE

    INTEGER meqn,mbc,mx,my,maux
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

    INTEGER i,j, mq
    DOUBLE PRECISION xi,yj

    DO mq = 1,meqn
        DO  i = 1-mbc,mx+mbc
            xi = xlower + (i-0.5d0)*dx
            DO j = 1-mbc,my+mbc
                yj = ylower + (j-0.5d0)*dy
                IF (xi .LT. 0.5d0) THEN
                    q(mq,i,j) = 1.d0
                ELSE
                    q(mq,i,j) = 0.d0
                ENDIF
            ENDDO
        ENDDO
    ENDDO


END SUBROUTINE clawpack5_qinit
