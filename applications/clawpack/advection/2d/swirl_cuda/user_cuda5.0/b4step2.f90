SUBROUTINE cudaclaw5_b4step2(mbc,mx,my,meqn,q, &
  xlower,ylower,dx,dy,time,dt,maux,aux)
  IMPLICIT NONE

  INTEGER mbc,mx,my,maux,meqn
  DOUBLE PRECISION xlower,ylower,dx,dy,time,dt

  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION tperiod, pi2
  COMMON /comvt/ tperiod,pi2

  INTEGER i,j
  DOUBLE PRECISION xll,yll, psi, vt

  IF (tperiod .EQ. 0.d0) THEN
     !! # special case --- indication that velocities specified in
     !! # setaux should be used for all time.
     RETURN
  ENDIF

  vt = dcos(pi2*(time+dt/2.d0)/tperiod)

  DO i=1-mbc,mx+mbc  
     DO  j=1-mbc,my+mbc
        !! # coordinates of lower left corner of grid cell:
        xll = xlower + (i-1)*dx
        yll = ylower + (j-1)*dy

        !! # difference stream function psi to get normal velocities:
        aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
        aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx

        !! # multiply by time-factor:
        aux(1,i,j) = vt * aux(1,i,j)
        aux(2,i,j) = vt * aux(2,i,j)  
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE cudaclaw5_b4step2
