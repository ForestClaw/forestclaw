SUBROUTINE filament_clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, my, maux
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER i, j, mq, blockno, fc2d_clawpack5_get_block
  DOUBLE PRECISION xlow, ylow, w

  blockno = fc2d_clawpack5_get_block()

  DO mq = 1,meqn
     DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
           xlow = xlower + (i-1)*dx
           ylow = ylower + (j-1)*dy
           CALL cellave2(blockno,xlow,ylow,dx,dy,w)
           q(mq,i,j) = w
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE filament_clawpack5_qinit
