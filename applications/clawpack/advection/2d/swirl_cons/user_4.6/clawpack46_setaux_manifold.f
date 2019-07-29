      subroutine clawpack46_setaux_manifold(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,blockno,
     &      xp,yp,zp,
     &      area, edgelengths,xnormals,ynormals,surfnormals)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision    xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision    ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

      integer i,j, k
      double precision dxdy, t, u, v, nv(2)

c     # Set velocity field;  For a manifold, we should project out
c     # the surface normals.

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # coordinates of lower left corner of grid cell:
            call velocity_field(xp(i,j),yp(i,j),u,v)
            aux(i,j,1) = u
            aux(i,j,2) = v
         enddo
      enddo

      dxdy = dx*dy

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,3) = area(i,j)/dxdy
         enddo
      enddo

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,4) = edgelengths(i,j,1)/dy
         enddo
      enddo

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,5) = edgelengths(i,j,2)/dx
         enddo
      enddo

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            do k = 1,2
               aux(i,j,5+k) = xnormals(i,j,k)
               aux(i,j,7+k) = ynormals(i,j,k)
            enddo
         enddo
      enddo

      return
      end

