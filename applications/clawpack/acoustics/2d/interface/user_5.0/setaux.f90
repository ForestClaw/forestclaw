SUBROUTINE clawpack5_setaux(mbc,mx,my, &
     xlower,ylower,dx,dy,maux,aux)
  IMPLICIT NONE

  INTEGER mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy

  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION rhol, cl, rhor, cr
  COMMON /comaux/ rhol,cl,rhor,cr

  INTEGER i,j
  DOUBLE PRECISION xcell

  DO j=1-mbc,my+mbc
     DO i=1-mbc,mx+mbc
        xcell = xlower + (i-0.5d0)*dx
        IF (xcell .LT. 0.5d0) THEN
           aux(1,i,j) = rhol
           aux(2,i,j) = cl
        ELSE
           aux(1,i,j) = rhor
           aux(2,i,j) = cr
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_setaux
