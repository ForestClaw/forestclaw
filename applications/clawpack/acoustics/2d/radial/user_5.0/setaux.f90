SUBROUTINE clawpack5_setaux_manifold(mbc,mx,my, &
        xlower,ylower,dx,dy,maux,aux, &
        xnormals,ynormals,edgelengths,area)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION        area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
  DOUBLE PRECISION     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION  edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)

  DOUBLE PRECISION rho,bulk,cc,zz
  COMMON /cparam/ rho,bulk,cc,zz

  INTEGER i,j

  DO  j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        aux(1,i,j) = xnormals(i,j,1)
        aux(2,i,j) = xnormals(i,j,2)
        aux(3,i,j) = edgelengths(i,j,1)/dy
        aux(4,i,j) = ynormals(i,j,1)
        aux(5,i,j) = ynormals(i,j,2)
        aux(6,i,j) = edgelengths(i,j,2)/dx
        aux(7,i,j) = area(i,j)/(dx*dy)
        aux(8,i,j) = cc
        aux(9,i,j) = zz
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_setaux_manifold
