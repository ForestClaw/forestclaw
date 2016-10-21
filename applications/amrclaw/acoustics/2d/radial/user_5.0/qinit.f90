SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux)
  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DIMENSION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

  pi = 4.d0*datan(1.d0)
  width = 0.2d0

  DO  i=1-mbc,mx+mbc
     xcell = xlower + (i-0.5d0)*dx
     DO j=1-mbc,my+mbc
        ycell = ylower + (j-0.5d0)*dy
        r = dsqrt(xcell**2 + ycell**2)

        IF (dabs(r-0.5d0) .LE. width) THEN
           pressure = 1.d0 + dcos(pi*(r - 0.5d0)/width)
        ELSE
           pressure = 0.d0
        ENDIF
        q(1,i,j) = pressure
        q(2,i,j) = 0.d0
        q(3,i,j) = 0.d0
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack5_qinit
