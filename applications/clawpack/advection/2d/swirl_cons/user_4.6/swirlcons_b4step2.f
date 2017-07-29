      subroutine swirlcons_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j, example
      double precision tperiod, pi2, vt, xc,yc, ucc, vcc, xe,ye

      common /comvt/ tperiod,pi2
      common /comex/ example
c
      if (tperiod .eq. 0.d0) then
          vt = 1.d0
      else
          vt = cos(pi2*(time+dt/2.d0)/tperiod)
      endif

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            xe = xlower + (i-1)*dx
            ye = ylower + (j-1)*dy

c           # Cell-centered velocity
            aux(i,j,1) = vt*ucc(xc,yc)
            aux(i,j,2) = vt*vcc(xc,yc)
         enddo
      enddo

      return
      end

      double precision function ucc(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

c      # 2d : Swirl example (computed from streamfunction)           
      ucc = 2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp))

      ucc = cos(pi*xp)**2 - 0.5

      return
      end

      double precision function vcc(xp,yp)
      implicit none
 
      double precision xp,yp,pi

      common /compi/ pi


c     # 2d examples
      vcc = -2*sin(pi*yp)**2 * sin(pi*xp) * cos(pi*xp)

c     # For a non-divergence free example
      vcc = vcc + 0.2d0*sin(2*pi*yp)*sin(2*pi*xp)

      vcc = sin(pi*yp)**2 - 0.5


      return
      end
