SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower, &
                           dx,dy,dz,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,mz,maux
  DOUBLE PRECISION xlower,ylower,zlower,dx,dy,dz
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)

  INTEGER i,j,k
  DOUBLE PRECISION xi,yj

  DO k = 1-mbc,mz+mbc
    DO  i = 1-mbc,mx+mbc
       xi = xlower + (i-0.5d0)*dx
       DO j = 1-mbc,my+mbc
          yj = ylower + (j-0.5d0)*dy
          IF (xi .LT. 0.5d0) THEN
             q(1,i,j,k) = 1.d0
          ELSE
             q(1,i,j,k) = 0.d0
          ENDIF
       ENDDO
    ENDDO
  ENDDO

END SUBROUTINE clawpack5_qinit
