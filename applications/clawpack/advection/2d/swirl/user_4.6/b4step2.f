      subroutine clawpack46_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision vt, xll,yll, psi, pij, vtx,vty

      double precision tperiod, pi2
      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

      vt = cos(pi2*(time+dt/2.d0)/tperiod)
      vty = vt/dy
      vtx = vt/dx

      do j = 1-mbc,my+mbc
         yll = ylower + (j-1)*dy
         do i = 1-mbc,mx+mbc
             xll = xlower + (i-1)*dx
             
             pij = psi(xll,yll)
             aux(i,j,1) = vty*(psi(xll, yll+dy) - pij)
             aux(i,j,2) = -vtx*(psi(xll+dx, yll) - pij)
        enddo
      enddo

      return
      end
