SUBROUTINE clawpack5_b4step2(mbc,mx,my,meqn,q, &
  xlower,ylower,dx,dy,time,dt,maux,aux)
  IMPLICIT NONE

  INTEGER :: mbc,mx,my,maux,meqn
  DOUBLE PRECISION :: xlower,ylower,dx,dy,time,dt

  DOUBLE PRECISION :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION :: pi, pi2
  COMMON /compi/ pi, pi2

  DOUBLE PRECISION :: tperiod
  COMMON /comvt/ tperiod

  INTEGER :: i,j
  DOUBLE PRECISION :: xll,yll, psi, vt, vtx,vty,pij

  IF (tperiod .EQ. 0.d0) THEN
     !! # special case --- indication that velocities specified in
     !! # setaux should be used for all time.
     RETURN
  ENDIF

  vt = dcos(pi2*(time+dt/2.d0)/tperiod)
  vtx = vt/dx
  vty = vt/dy

  DO  j=1-mbc,my+mbc
    DO i=1-mbc,mx+mbc
        !! # coordinates of lower left corner of grid cell:
        xll = xlower + (i-1)*dx
        yll = ylower + (j-1)*dy

        !! # difference stream function psi to get normal velocities:
        pij = psi(xll,yll)
        aux(1,i,j) =  vty*(psi(xll, yll+dy) - pij)
        aux(2,i,j) = -vtx*(psi(xll+dx, yll) - pij)

     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack5_b4step2
