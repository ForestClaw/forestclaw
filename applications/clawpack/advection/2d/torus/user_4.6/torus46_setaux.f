      subroutine torus46_setaux(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,blockno,
     &      area, edgelengths,
     &      xnormals,ynormals,surfnormals)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, k
      double precision dxdy, xc, yc, nv(3), t, u, v, w
      double precision u1, u2, v1, v2, e1,e2,e3,e4,udiv

      include "metric_terms.i"

c     # ----------------------------------------------------------------
c     # 1        capacity
c     # 2-3      Edge velocities
c     # 4-6      Cell-centered velocities (x,y,z)
c     # 7-8      Edge lengths (x-face, y-face)
c     # 9-11     x-face normals
c     # 12-14    y-face normals
c     # ----------------------------------------------------------------

      t = 0
      dxdy = dx*dy

c     # Capacity : entry (1)
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

c     # Edge velocities : entries (2,3)      
      call torus46_set_edge_velocities(mx,my,mbc,dx,dy,
     &       blockno,xlower,ylower,aux,maux)

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            do k = 1,3
               nv(k) = surfnormals(i,j,k)
            enddo

            call torus_center_velocity(blockno,xc,yc,t, nv, 
     &                                  u, v, w)
            aux(i,j,4) = u
            aux(i,j,5) = v
            aux(i,j,6) = w
         enddo
      enddo
100   format(4F16.8,E16.8)      

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face and y-face edge lengths (6,7)      
            aux(i,j,7) = edgelengths(i,j,1)/dy
            aux(i,j,8) = edgelengths(i,j,2)/dx

c           # Normals (9,10,11) and (12,13,14)
            do k = 1,3
               aux(i,j,9  + k-1) = xnormals(i,j,k)
               aux(i,j,12 + k-1) = ynormals(i,j,k)
            enddo

         enddo
      enddo

      return
      end


      subroutine torus46_set_edge_velocities(mx,my,mbc,
     &      dx,dy,blockno,xlower,ylower,aux,maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower

      double precision xd1(2),xd2(2), t
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      t = 0.d0  !! not used
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xlower + (i-1)*dx
            xd1(2) = ylower + j*dy

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy

            call torus_edge_velocity(blockno,xd1,xd2,dy,vn,t)
            aux(i,j,2) = vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # y-faces
            xd1(1) = xlower + i*dx
            xd1(2) = ylower + (j-1)*dy

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy

            call torus_edge_velocity(blockno,xd1,xd2,dx,vn,t)
            aux(i,j,3) = -vn
         enddo
      enddo

      end

