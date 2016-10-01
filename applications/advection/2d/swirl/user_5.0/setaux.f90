SUBROUTINE clawpack5_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
  IMPLICIT NONE

  INTEGER mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION psi

  INTEGER i,j
  DOUBLE PRECISION xll, yll
  !! # constant velocities which are used if tperiod=0 is specified
  !! # in setprob.data

  DO  i=1-mbc,mx+mbc
     DO  j=1-mbc,my+mbc

        !!  # coordinates of lower left corner of grid cell:
        xll = xlower + (i-1)*dx
        yll = ylower + (j-1)*dy

        !!  # difference stream function psi to get normal velocities:
        aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
        aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack5_setaux
