      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i, j
      double precision xc,yc,u,v
      integer blockno, fc2d_clawpack46_get_block

      blockno = fc2d_clawpack46_get_block()

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            call velocity_field(xc,yc,u,v)
            aux(i,j,1) = u
            aux(i,j,2) = v
 
         enddo
      enddo

      return
      end



