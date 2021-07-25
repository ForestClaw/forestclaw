SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my, &
     xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy

  DOUBLE PRECISION    q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER blockno, fc2d_clawpack5_get_block

  DOUBLE PRECISION hin,hout
  COMMON /comic/ hin,hout

  INTEGER i,j
  DOUBLE PRECISION xc,yc

  blockno = fc2d_clawpack5_get_block()

  DO i = 1-mbc,mx+mbc
     xc = xlower + (i-0.5)*dx
     DO j = 1-mbc,my+mbc
        yc = ylower + (j-0.5)*dy
        q(1,i,j) = 0.1d0 + exp(-200.d0*(xc**2 + yc**2))
        q(2,i,j) = 0.d0
        q(3,i,j) = 0.d0
     ENDDO
  ENDDO
END SUBROUTINE clawpack5_qinit
