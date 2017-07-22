      subroutine swirlcons_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j, example
      double precision tperiod, pi2, vt, xc,yc, ucc, vcc

      common /comvt/ tperiod,pi2
      common /comex/ example
c
      if (example .le. 3) then
          vt = 1.d0
      else
          if (tperiod .eq. 0.d0) then
              vt = 1.d0
          else
              vt = cos(pi2*(time+dt/2.d0)/tperiod)
          endif
      endif

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # Cell-centered velocity
            aux(i,j,1) = vt*ucc(xc,yc)
            aux(i,j,2) = vt*vcc(xc,yc)
         enddo
      enddo

      return
      end

      double precision function ucc(xp,yp)
      implicit none

      double precision xp,yp,pi, c,rp2,ul,ur,uavg
      integer example
      common /compi/ pi
      common /comex/ example

      if (example .eq. 1) then
c         # 1d : Piecewise constant velocity (-0.5, 0.5)
          if (xp .le. 0.5d0) then
              ucc = -0.5d0
          else
              ucc = 0.5d0
          endif
      elseif (example .eq. 2) then
c         # 1d : u = cos(2*pi*xp)       
          ucc = cos(2*pi*xp) + 1.d0            
      elseif (example .eq. 3) then
c         # 1d : More complicated 1d example
          ucc = 0.1*sin(2*pi*xp)*sin(16*pi*xp)
      elseif ((example .eq. 4) .or. (example .eq. 5)) then
c         # 2d : Swirl example (computed from streamfunction)           
          ucc = 2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp))
          if (example .eq. 5) then
c             # Add non-zero divergence            
              rp2 = (xp-0.5d0)**2 + (yp-0.5d0)**2
              rp2 = (xp-0.5d0)**2
              c = exp(-350.d0*rp2)
              ucc = ucc + 5*c
          endif
      else
          write(6,*) 'b4step2 : example not set'
          stop
      endif

      return
      end

      double precision function vcc(xp,yp)
      implicit none

      double precision xp,yp,pi, rp2, c
      integer example

      common /compi/ pi
      common /comex/ example


      if (example .le. 3) then
c         # 1d examples
          vcc = 0        
      elseif ((example .eq. 4) .or. (example .eq. 5)) then
c         # 2d examples
          vcc = -2*((sin(pi*yp))**2 * sin(pi*xp) * cos(pi*xp))
          if (example .eq. 5) then
c             # Add non-zero divergence              
              rp2 = (xp-0.5d0)**2 + (yp-0.5d0)**2
              rp2 = (xp-0.5d0)**2
              c = exp(-350.d0*rp2)
              vcc = vcc + 5*c
          endif
      endif

      return
      end
