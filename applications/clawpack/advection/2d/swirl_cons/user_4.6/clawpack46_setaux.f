      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i, j
      double precision xc, yc, pi,s
      integer example

      common /compi/ pi
      common /comex/ example

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            s = sqrt(2.d0)
            if (example .eq. 1) then
               aux(i,j,1) = s*(cos(pi*xc)**2 + 0.5d0)
               aux(i,j,2) = s*(sin(pi*yc)**2 + 0.5d0)
            else
               aux(i,j,1) = s*(cos(pi*xc)**2 - 0.5d0)
               aux(i,j,2) = s*(sin(pi*yc)**2 - 0.5d0)
            endif

         enddo
      enddo

      return
      end
