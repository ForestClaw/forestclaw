      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)

      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, maux)

      double precision xe,xc, pi, dx1, dx2, a
      integer i, example

      common /compi/ pi
      common /comex/ example

c     # aux(i,1) : edge velocity (u_{i-1/2})
c     # aux(i,2) : cell-centered derivative of velocity (for src2)

      do i = 1-mbc,mx+mbc
          xe = xlower + (i-1)*dx
          xc = xlower + (i-0.5)*dx

          if (example .eq. 1) then
c             # Velocity should have a delta function at x=0.5
              if (x .le. 0.5d0) then
                  aux(i,1) = -0.5d0
                  aux(i,2) = 0
              else
                  aux(i,1) = 0.5d0
                  aux(i,2) = 0            
              endif
          elseif (example .eq. 2) then
              aux(i,1) = cos(2*pi*xe)
              aux(i,2) = -2*pi*sin(2*pi*xc)
          elseif (example .eq. 3) then
              aux(i,1) = 0.1d0*sin(2*pi*xe)*sin(16*pi*xe)
              dx1 =  2*pi*cos(2*pi*xc)*sin(16*pi*xc)
              dx2 = 16*pi*sin(2*pi*xc)*cos(16*pi*xc)
              aux(i,2) = 0.1*(dx1 + dx2)
          elseif (example .eq. 4) then
              a = 0.01d0
              aux(i,1) = tanh((xe-0.5d0)/a)
              aux(i,2) = (1-tanh((xc-0.5d0)/a)**2)/a
          elseif (example .eq. 5) then
              a = 0.1d0
              aux(i,1) = -tanh((xe-0.5d0)/a)
              aux(i,2) = -(1-tanh((xc-0.5d0)/a)**2)/a
          endif
      enddo

      return
      end
