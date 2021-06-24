SUBROUTINE user5_setaux_manifold(mbc,mx,my, &
        xlower,ylower,dx,dy,maux,aux, &
        xnormals,xtangents,ynormals,ytangents, &
        surfnormals,area)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  !! Still hoping to be able to include this using 'fclaw2d_metric_terms.i'
  DOUBLE PRECISION         area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
  DOUBLE PRECISION  surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
  DOUBLE PRECISION     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION    xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
  DOUBLE PRECISION    ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)


  INTEGER i,j

  DO  j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        aux(1,i,j) = area(i,j)/(dx*dy)

        aux(2,i,j) = xnormals(i,j,1)
        aux(3,i,j) = xnormals(i,j,2)
        aux(4,i,j) = xnormals(i,j,3)

        aux(5,i,j) = xtangents(i,j,1)
        aux(6,i,j) = xtangents(i,j,2)
        aux(7,i,j) = xtangents(i,j,3)

        aux(8,i,j)  = ynormals(i,j,1)
        aux(9,i,j)  = ynormals(i,j,2)
        aux(10,i,j) = ynormals(i,j,3)
        
        aux(11,i,j) = ytangents(i,j,1)
        aux(12,i,j) = ytangents(i,j,2)
        aux(13,i,j) = ytangents(i,j,3)

        aux(14,i,j) = surfnormals(i,j,1)
        aux(15,i,j) = surfnormals(i,j,2)
        aux(16,i,j) = surfnormals(i,j,3)
     ENDDO
  ENDDO

END SUBROUTINE user5_setaux_manifold
