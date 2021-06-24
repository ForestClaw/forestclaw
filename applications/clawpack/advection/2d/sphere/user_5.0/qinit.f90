SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, my, maux
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION xlow, ylow, w
  INTEGER blockno, fc2d_clawpack5_get_block

  blockno = fc2d_clawpack5_get_block()

  DO j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        xlow = xlower + (i-1)*dx
        ylow = ylower + (j-1)*dy
        CALL cellave2(blockno,xlow,ylow,dx,dy,w)
        q(1,i,j) = w
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_qinit
