      subroutine clawpack46_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision tperiod, pi2, vt, xll,yll, psi, u, v

      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

      vt = cos(pi2*(time+dt/2.d0)/tperiod)

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
C             aux(i,j,1) = (psi(xll, yll+dy) - psi(xll,yll)) / dy
C             aux(i,j,2) =  -(psi(xll+dx, yll) - psi(xll,yll)) / dx
            
c           # Cell-centered velocity
            aux(i,j,1) = u(xll,yll)
            aux(i,j,2) = v(xll,yll)

c           # multiply by time-factor:
            aux(i,j,1) = vt * aux(i,j,1)
            aux(i,j,2) = vt * aux(i,j,2)
C             aux(i,j,1) = 1.0
C             aux(i,j,2) = 0.0
         enddo
      enddo

      return
      end

      double precision function u(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

      u = 2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp))

      return
      end

      double precision function v(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

      v = -2*((sin(pi*yp))**2 * sin(pi*xp) * cos(pi*xp))

      return
      end
