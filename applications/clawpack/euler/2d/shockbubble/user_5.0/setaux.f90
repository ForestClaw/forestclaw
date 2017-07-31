SUBROUTINE clawpack5_setaux(mbc,mx,my, &
     xlower,ylower,dx,dy,maux,aux)
  IMPLICIT NONE

  INTEGER mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION yj

  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        yj = ylower + (j-0.5d0)*dy
        aux(1,i,j) = yj
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_setaux
