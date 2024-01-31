SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, my, maux
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION xc,yc,xp,yp,zp
  DOUBLE PRECISION gaussian_sum

  INTEGER*8 cont, fclaw_map_get_context
  INTEGER blockno, fc2d_clawpack5_get_block

  cont = fclaw_map_get_context()
  blockno = fc2d_clawpack5_get_block()

  DO j = 1-mbc,my+mbc
        yc = ylower + (j-0.5d0)*dy
     DO i = 1-mbc,mx+mbc
        xc = xlower + (i-0.5d0)*dx
        CALL fclaw_map_2d_c2m(cont, &
             blockno,xc,yc,xp,yp,zp)
        q(1,i,j) = gaussian_sum(xp,yp,zp)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_qinit
