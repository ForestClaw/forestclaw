SUBROUTINE torus5_compute_error(blockno, mx,my,mbc,meqn, &
     dx,dy,xlower,ylower,t, q,error)
  IMPLICIT NONE

  INTEGER mx,my,mbc,meqn, blockno
  DOUBLE PRECISION dx, dy, xlower, ylower, t
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j,m
  DOUBLE PRECISION xc,yc, qexact

  !! # Assume a single field variable only
  DO j = 1,my
     yc = ylower + (j-0.5)*dy
     DO i = 1,mx
        xc = xlower + (i-0.5)*dx
        error(1,i,j) = q(1,i,j) - qexact(blockno,xc,yc,t);
     ENDDO
  ENDDO


END SUBROUTINE torus5_compute_error
