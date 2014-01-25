      subroutine b4step2(mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i, j
      double precision tperiod, pi2, vt, xll,yll, psi

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
            aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
c
c           # multiply by time-factor:
            aux(1,i,j) = vt * aux(1,i,j)
            aux(2,i,j) = vt * aux(2,i,j)
         enddo
      enddo

      return
      end
