      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Shallow water with radial dam break problem, h = hin inside
c     # circle specified in fdisc.f
c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision dx, dy, xlower, ylower
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision xc, yc

      do i = 1-mbc,mx+mbc
          xc = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
              yc = ylower + (j-0.5d0)*dy
              q(i,j,1) = 0.1d0 + exp(-200.d0*(xc**2 + yc**2))
              q(i,j,2) = 0.d0
              q(i,j,3) = 0.d0
          end do
      end do

      return
      end
