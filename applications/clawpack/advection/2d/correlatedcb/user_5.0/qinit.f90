SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, my, maux
  DOUBLE PRECISION xlower, ylower, dx, dy
  INTEGER*8 cont, get_context
  INTEGER blockno, fc2d_clawpack5_get_block
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION xc,yc,xp,yp,zp
  DOUBLE PRECISION cosine_bell_sum

  DOUBLE PRECISION a,b

  cont = get_context()
  blockno = fc2d_clawpack5_get_block()

  a = -0.8d0
  b = 0.9d0

  DO j = 1-mbc,my+mbc
        yc = ylower + (j-0.5d0)*dy
     DO i = 1-mbc,mx+mbc
        xc = xlower + (i-0.5d0)*dx
        CALL fclaw2d_map_c2m(cont, &
             blockno,xc,yc,xp,yp,zp)
        q(1,i,j) = cosine_bell_sum(xp,yp,zp)
        IF (meqn .EQ. 2) THEN
           !! # Set non-linear relationship between two tracers
           q(2,i,j) = a*q(1,i,j)**2 + b
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_qinit
