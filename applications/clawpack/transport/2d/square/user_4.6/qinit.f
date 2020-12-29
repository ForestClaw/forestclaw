      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision xc,yc,dsum, rc, delta

      do i = 1-mbc,mx+mbc
          xc = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
              yc = ylower + (j-0.5d0)*dy
c              rc = sqrt((xc-0.5d0)**2 + (yc-0.5d0)**2)
c              q(i,j,1) = delta(rc)
               q(i,j,1) = 0
          enddo
      enddo

      return
      end


      double precision function delta2(x)
      implicit none

      double precision x

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision t0

      t0 = 1e-4
      delta2 = exp(-x**2/(4*t0))/(4*pi*t0)

      return

      end function delta2


      
