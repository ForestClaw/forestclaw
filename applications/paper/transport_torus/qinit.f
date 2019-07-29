       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      double precision pi, pi2
      common /compi/ pi, pi2

      integer i,j
      double precision xc,yc, th, gm

      do j = 1-mbc,my+mbc
          yc = ylower + (j-0.5)*dy
          gm = pi2*yc
          do i = 1-mbc,mx+mbc     
              xc = xlower + (i-0.5)*dx
              th = pi2*xc
              if (abs(th-pi) .le. pi/2 .and. abs(gm-pi) .le. pi/2) then
                  q(1,i,j) = 1.d0
              else
                  q(1,i,j) = 0
              endif
          enddo
      enddo

      return
      end
