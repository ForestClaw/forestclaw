      subroutine user46_setaux_manifold(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,
     &      xnormals,ynormals,edgelengths,area)
      implicit none

      integer mx,my,mbc,maux
      double precision xlower,ylower,dx,dy
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision rho,bulk,cc,zz
      COMMON /cparam/ rho,bulk,cc,zz

      integer i,j

      include "metric_terms.i"

      do  j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            aux(i,j,1) = xnormals(i,j,1)
            aux(i,j,2) = ynormals(i,j,2)
            aux(i,j,3) = edgelengths(i,j,1)/dy
            aux(i,j,4) = ynormals(i,j,1)
            aux(i,j,5) = ynormals(i,j,2)
            aux(i,j,6) = edgelengths(i,j,2)/dx
            aux(i,j,7) = area(i,j)/(dx*dy)
            aux(i,j,8) = cc
            aux(i,j,9) = zz
         enddo
      enddo

      return
      end
