      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i, j
      double precision xc, yc, ucc, pi

      common /compi/ pi

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           ucc = 2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp))
c           ucc = cos(2*pi*xp)

            ucc = 0.1*sin(2*pi*xc)*sin(16*pi*xc)

            aux(i,j,1) = ucc
            aux(i,j,2) = 0.d0
         enddo
      enddo

      return
      end
