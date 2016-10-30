SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my, &
     xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower, dx,dy
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION xi,yj,r, p0

  DO i=1-mbc,mx+mbc
     xi = xlower + (i-0.5d0)*dx
     DO j=1-mbc,my+mbc
        yj = ylower + (j-0.5d0)*dy
        r = dsqrt((xi-0.25d0)**2 + (yj-0.4d0)**2)
        q(1,i,j) = p0(r)
        q(2,i,j) = 0.d0
        q(3,i,j) = 0.d0
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack5_qinit
