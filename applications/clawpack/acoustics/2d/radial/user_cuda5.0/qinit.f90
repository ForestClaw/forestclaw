SUBROUTINE cudaclaw5_qinit(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  DOUBLE PRECISION xc,yc,xp,yp,zp, width, r, pressure

  INTEGER i,j

  INTEGER blockno, fc2d_clawpack5_get_block
  INTEGER*8 cont, get_context
  LOGICAL fclaw2d_map_is_used

  DOUBLE PRECISION pi
  COMMON /compi/ pi

  blockno = fc2d_clawpack5_get_block()
  cont = get_context()

  width = 0.2d0

  DO  i = 1-mbc,mx+mbc
     xc = xlower + (i-0.5d0)*dx
     DO j = 1-mbc,my+mbc
        yc = ylower + (j-0.5d0)*dy
        IF (fclaw2d_map_is_used(cont)) THEN
           CALL fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
           r = SQRT(xp**2 + yp**2)
        ELSE
           r = SQRT(xc**2 + yc**2)
        ENDIF
        IF (ABS(r-0.5d0) .LE. width) THEN
           pressure = 1.d0 + COS(pi*(r - 0.5d0)/width)
        ELSE
           pressure = 0.d0
        ENDIF
        q(1,i,j) = pressure
        q(2,i,j) = 0.d0
        q(3,i,j) = 0.d0

     ENDDO
  ENDDO
END SUBROUTINE cudaclaw5_qinit
