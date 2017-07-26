      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)

      implicit none
      integer maxmx, mbc,mx,maux
      double precision dx, xlower
      double precision aux(1-mbc:maxmx+mbc, maux)

      double precision xc, pi, a, y1, y2, a1, a2
      integer i, example

      common /compi/ pi
      common /comex/ example

c     # aux(i,1) : edge velocity (u_{i-1/2})
c     # aux(i,2) : cell-centered derivative of velocity (for src2)

      do i = 1-mbc,mx+mbc
          xc = xlower + (i-0.5)*dx

          if (example .eq. 1) then
c             # Velocity should have a delta function at x=0.5
              if (xc .le. 0.5d0) then
                  aux(i,1) = -0.5d0
                  aux(i,2) = 0
              else
                  aux(i,1) = 0.5d0
                  aux(i,2) = 0            
              endif
          elseif (example .eq. 2) then
              aux(i,1) = cos(2*pi*xc) + 2.d0
          elseif (example .eq. 3) then
              a = 0.1d0
              aux(i,1) = a*sin(2*pi*xc)*sin(16*pi*xc)
          elseif (example .eq. 4) then
              a1 = 0.025d0
              a2 = 0.04d0
              y1 = 0.5*tanh((xc-0.65d0)/a1) + 1
              y2 = 0.5*tanh((0.35d0-xc)/a2) + 1
              aux(i,1) = y1 + y2 - 1.5
          elseif (example .eq. 5) then
              a = 0.01d0
              aux(i,1) = tanh((xc-0.5d0)/a)
          elseif (example .eq. 6) then
              a = 0.1d0
              aux(i,1) = -tanh((xc-0.5d0)/a)
          endif
      enddo

      return
      end
