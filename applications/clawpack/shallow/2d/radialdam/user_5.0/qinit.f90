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
  DOUBLE PRECISION xlow,ylow, win

  blockno = fc2d_clawpack5_get_block()

  DO i = 1-mbc,mx+mbc
     xlow = xlower + (i-1.d0)*dx
     DO j = 1-mbc,my+mbc
        ylow = ylower + (j-1.d0)*dy
        CALL cellave2(blockno,xlow,ylow,dx,dy,win)
        q(1,i,j) = hin*win + hout*(1.d0-win)
        q(2,i,j) = 0.d0
        q(3,i,j) = 0.d0
     ENDDO
  ENDDO
END SUBROUTINE clawpack5_qinit
