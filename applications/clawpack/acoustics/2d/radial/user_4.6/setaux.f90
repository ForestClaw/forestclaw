SUBROUTINE clawpack46_setaux_manifold(mbc,mx,my, &
        xlower,ylower,dx,dy,maux,aux, &
        xnormals,ynormals,edgelengths,area)
  IMPLICIT NONE

  INTEGER :: mx,my,mbc,maux
  DOUBLE PRECISION :: xlower,ylower,dx,dy
  DOUBLE PRECISION :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

  !! Still hoping to be able to include this using 'fclaw2d_metric_terms.i'
  DOUBLE PRECISION ::         area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
  DOUBLE PRECISION ::    xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION ::    ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION :: edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)

  DOUBLE PRECISION rho,bulk,cc,zz
  COMMON /cparam/ rho,bulk,cc,zz

  INTEGER i,j

  DO  j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        aux(i,j,1) = xnormals(i,j,1)
        aux(i,j,2) = xnormals(i,j,2)
        aux(i,j,3) = edgelengths(i,j,1)/dy
        aux(i,j,4) = ynormals(i,j,1)
        aux(i,j,5) = ynormals(i,j,2)
        aux(i,j,6) = edgelengths(i,j,2)/dx
        aux(i,j,7) = area(i,j)/(dx*dy)
        aux(i,j,8) = cc
        aux(i,j,9) = zz
     ENDDO
  ENDDO

END SUBROUTINE clawpack46_setaux_manifold
