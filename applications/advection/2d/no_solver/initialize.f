      subroutine initialize(mx,my,meqn,mbc,xlower,ylower,dx,dy,q)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

      integer i, j
      double precision xlow, ylow,wl

      do i = 1-mbc,mx+mbc
         xlow = xlower + (i-1)*dx
         do j = 1-mbc,my+mbc
            ylow = ylower + (j-1)*dy

            call cellave2(xlow,ylow,dx,dy,wl)

            q(i,j,1) = wl
            if (meqn > 1) then
               q(i,j,2) = 1-q(i,j,1)
            endif
         enddo
      enddo

      return
      end

      double precision function fdisc(x,y)
      implicit none

      double precision x,y, r

      r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

      fdisc = r-0.25d0

      end
