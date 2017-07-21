      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)

      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, *)

      double precision xe,xc, pi, dx1, dx2
      integer i

      common /compi/ pi


      do i = 1-mbc,mx+mbc
c          x = xlower + (i-0.5)*dx
          xe = xlower + (i-1)*dx
          xc = xlower + (i-0.5)*dx

          aux(i,1) = cos(2*pi*xe)
          aux(i,2) = -2*pi*sin(2*pi*xc)

c          aux(i,1) = 0.1d0*sin(2*pi*xe)*sin(16*pi*xe)
c          dx1 =  2*pi*cos(2*pi*xc)*sin(16*pi*xc)
c          dx2 = 16*pi*sin(2*pi*xc)*cos(16*pi*xc)
c          aux(i,2) = 0.1*(dx1 + dx2)

c          if (x .le. 0.5d0) then
c              aux(i,1) = -0.5d0
c          else
c              aux(i,1) = 0.5d0
c          endif
      enddo

c
      return
      end
