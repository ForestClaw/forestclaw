      subroutine torus46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)

      implicit none

      integer maxmx, maxmy, mbc, mx,my, maux
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer blockno, fc2d_clawpack46_get_block

      blockno = fc2d_clawpack46_get_block()

c     # Edge velocities : entries (2-3)      
      call torus46_set_edge_velocities(mx,my,mbc,dx,dy,
     &             blockno,xlower,ylower,aux,maux)

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

      integer*8 cont, fclaw_map_get_context

      integer i,j
      double precision vn

      cont = fclaw_map_get_context()

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face - upper vertex
            xc = xlower + (i-1)*dx
            yc = ylower + j*dy

c           # Map the brick to a unit square      
            call fclaw_map_2d_brick2c(cont, blockno,xc,yc,xc1,yc1,zc1)

c           # x-face - lower vertex
            xc = xlower + (i-1)*dx
            yc = ylower + (j-1)*dy

c           # Map the brick to a unit square      
            call fclaw_map_2d_brick2c(cont, blockno,xc,yc,xc2,yc2,zc2)

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
            call fclaw_map_2d_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

c           # y-face - left vertex
            xc = xlower + (i-1)*dx
            yc = ylower + (j-1)*dy
            call fclaw_map_2d_brick2c(cont,blockno,xc,yc,xc2,yc2,zc2)

            call torus_edge_velocity(xc1,yc1,xc2,yc2,dx,vn)
            aux(i,j,3) = -vn
         enddo
      enddo

      end



