      subroutine swirlcons_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision tperiod, pi2, vt, xc,yc, psi, ucc, vcc

      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

      vt = cos(pi2*(time+dt/2.d0)/tperiod)
c      vt = 1.0

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # difference stream function psi to get normal velocities:
C             aux(i,j,1) = (psi(xll, yll+dy) - psi(xll,yll)) / dy
C             aux(i,j,2) =  -(psi(xll+dx, yll) - psi(xll,yll)) / dx
            
c           # Cell-centered velocity
            aux(i,j,1) = ucc(xc,yc)
            aux(i,j,2) = vcc(xc,yc)

c           # multiply by time-factor:
            aux(i,j,1) = vt * aux(i,j,1)
            aux(i,j,2) = vt * aux(i,j,2)
C             aux(i,j,1) = 1.0
C             aux(i,j,2) = 0.0
         enddo
      enddo

      return
      end

      double precision function ucc(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

c      ucc = 2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp))
      ucc = cos(2*pi*xp)

c      ucc = 0.1*sin(2*pi*xp)*sin(16*pi*xp)
c      if (xp .le. 0.5d0) then
c            ucc = -.5d0
c      else
c            ucc = .50
c      endif
c      if(xp .gt. 0) then
c            ucc = exp(-(xp-0.5)**2/(2*0.1**2))
c      else
c            ucc = 0
c      endif

C      if (xp .gt. 0.5) then
C            ucc = 0.5
C      else
C            ucc = 1.0
C      endif

c      if((xp .gt. 0.5 .and. xp .lt. 0.75)) then
c            ucc = 0.1
c      else if((xp .gt. 0.25 .and. xp .lt. 0.5)) then
c            ucc = -0.1
c      else
c            ucc = 0.0
c      endif

      return
      end

      double precision function vcc(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

c      vcc = -2*((sin(pi*yp))**2 * sin(pi*xp) * cos(pi*xp))
      vcc = 0

      return
      end
