SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, my, maux
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  DOUBLE PRECISION xlow, ylow, w, xc,yc, q0
  double precision xp,yp,zp, q0_physical
  INTEGER i,j

  INTEGER blockno, fc2d_clawpack46_get_block
  INTEGER*8 cont, get_context

  INTEGER initchoice
  COMMON /initchoice_comm/ initchoice

  blockno = fc2d_clawpack46_get_block()

  cont = get_context()


  DO j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        IF (initchoice .EQ. 0) THEN
           !! # Discontinuous solution
           xlow = xlower + (i-1)*dx
           ylow = ylower + (j-1)*dy
           CALL cellave2(blockno,xlow,ylow,dx,dy,w)
           q(1,i,j) = w
        ELSEIF (initchoice .EQ. 1) THEN
           !! # Smooth solution for computing the error
           xc = xlower + (i-0.5)*dx
           yc = ylower + (j-0.5)*dy
           call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
           q(1,i,j) = q0_physical(xp,yp,zp)
        ELSEIF (initchoice .EQ. 2) then
           q(1,i,j) = 1.d0
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_qinit
