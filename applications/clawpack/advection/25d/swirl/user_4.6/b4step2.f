      subroutine clawpack46_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j, k
      double precision tperiod, pi2, vt, xll,yll, psi

      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

      vt = cos(pi2*(time+dt/2.d0)/tperiod)

      do k = 1,2
         do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
c               # coordinates of lower left corner of grid cell:
                xll = xlower + (i-1)*dx
                yll = ylower + (j-1)*dy

c               # difference stream function psi to get normal velocities:
                if (k .eq. 1) then
                    aux(i,j,k) = (psi(xll, yll+dy) - psi(xll,yll)) / dy
                else
                    aux(i,j,k) = -(psi(xll+dx, yll) - psi(xll,yll)) / dx
                endif

c               # multiply by time-factor:
                aux(i,j,k) = vt * aux(i,j,k)
c                aux(i,j,2) = vt * aux(i,j,2)
           enddo
        enddo
      enddo

      return
      end
