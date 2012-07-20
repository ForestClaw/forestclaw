c     ============================================
      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,time,dt,maux,aux)
c     ============================================
c
c     # called before each call to step
c     # use to set time-dependent aux arrays or perform other tasks.
c
c     # make velocity time dependent, reversing flow.
c
c
      implicit none

      integer maxmx, maxmy, mbc, mx, my, meqn, maux
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision tperiod, pi2, vt, xll,yll, psi

      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

       vt = sin(pi2*(time+dt/2.d0)/tperiod)
c
       do i = 1-mbc,mx+mbc
          do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(i,j,1) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(i,j,2) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
c
c           # multiply by time-factor:
            aux(i,j,1) = vt * aux(i,j,1)
            aux(i,j,2) = vt * aux(i,j,2)
         enddo
      enddo

      return
      end
