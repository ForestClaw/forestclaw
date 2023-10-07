      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c
c     # Set initial conditions for q.
c     # Acoustics with smooth radially symmetric profile to test accuracy
c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx,dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision pi, width, xcell, ycell, r, pressure
      integer i,j

      pi = 4.d0*datan(1.d0)
      width = 0.2d0

      do i = 1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
             ycell = ylower + (j-0.5d0)*dy
             r = dsqrt(xcell**2 + ycell**2)

             if (dabs(r-0.5d0) .le. width) then
                 pressure = 1.d0 + dcos(pi*(r - 0.5d0)/width)
               else
                 pressure = 0.d0
             endif
             q(i,j,1) = pressure
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
          end do
      end do
       return
       end
