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
      double precision dxdy

      include "metric_terms.i"

c     # ----------------------------------------------------------------
c     # 1        capacity
c     # 2-3      Edge velocities
c     # 4-6      Cell-centered velocities (x,y,z)
c     # 7-8      Edge lengths (x-face, y-face)
c     # 9-11     x-face normals
c     # 12-14    y-face normals
c     # ----------------------------------------------------------------

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

      call torus46_set_center_velocities(mx,my,mbc,dx,dy,
     &       blockno,xlower,ylower,aux,maux, surfnormals)

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

      double precision xc,yc
      double precision xc1, yc1, zc1, xc2, yc2, zc2
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer*8 cont, get_context

      integer i,j
      double precision vn

      cont = get_context()

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face - upper vertex
            xc = xlower + (i-1)*dx
            yc = ylower + j*dy

c           # Map the brick to a unit square      
            call fclaw2d_map_brick2c(cont, blockno,xc,yc,xc1,yc1,zc1)

c           # x-face - lower vertex
            xc = xlower + (i-1)*dx
            yc = ylower + (j-1)*dy

c           # Map the brick to a unit square      
            call fclaw2d_map_brick2c(cont, blockno,xc,yc,xc2,yc2,zc2)

            call torus_edge_velocity(xc1,yc1,xc2,yc2,dy,vn)
            aux(i,j,2) = vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc

c           # Map (xc,yc,blockno) to [0,1]x[0,1]

c           # y-face - right vertex
            xc = xlower + i*dx
            yc = ylower + (j-1)*dy
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

c           # y-face - left vertex
            xc = xlower + (i-1)*dx
            yc = ylower + (j-1)*dy
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc2,yc2,zc2)

            call torus_edge_velocity(xc1,yc1,xc2,yc2,dx,vn)
            aux(i,j,3) = -vn
         enddo
      enddo

      end

      subroutine torus46_set_center_velocities(mx,my,mbc,
     &      dx,dy,blockno,xlower,ylower,aux,maux, surfnormals)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

      double precision xc,yc
      double precision xc1,yc1,zc1, nv(3), vel(3), vdotn, torus_dot


      integer*8 cont, get_context

      integer i,j, k

      cont = get_context()

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # This is not the torus mapping, but rather maps the brick to
c           # a unit square      
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            call torus_center_velocity(xc1,yc1,vel)

c           # subtract out normal components
            do k = 1,3
               nv(k) = surfnormals(i,j,k)
            enddo

            vdotn = torus_dot(vel,nv)

            aux(i,j,4) = vel(1) - vdotn*nv(1)
            aux(i,j,5) = vel(2) - vdotn*nv(2)
            aux(i,j,6) = vel(3) - vdotn*nv(3)
         enddo
      enddo

      end



