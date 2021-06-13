SUBROUTINE user46_setaux_manifold(mbc,mx,my, &
        xlower,ylower,dx,dy,maux,aux, &
        xnormals,xtangents,ynormals,ytangents, &
        surfnormals,area)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

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
        aux(i,j,1) = area(i,j)/(dx*dy)

        aux(i,j,2) = xnormals(i,j,1)
        aux(i,j,3) = xnormals(i,j,2)
        aux(i,j,4) = xnormals(i,j,3)
        aux(i,j,5) = xtangents(i,j,1)
        aux(i,j,6) = xtangents(i,j,2)
        aux(i,j,7) = xtangents(i,j,3)

        aux(i,j,8)  = ynormals(i,j,1)
        aux(i,j,9)  = ynormals(i,j,2)
        aux(i,j,10) = ynormals(i,j,3)
        aux(i,j,11) = ytangents(i,j,1)
        aux(i,j,12) = ytangents(i,j,2)
        aux(i,j,13) = ytangents(i,j,3)

        aux(i,j,14) = surfnormals(i,j,1)
        aux(i,j,15) = surfnormals(i,j,2)
        aux(i,j,16) = surfnormals(i,j,3)
     ENDDO
  ENDDO

END SUBROUTINE user46_setaux_manifold
