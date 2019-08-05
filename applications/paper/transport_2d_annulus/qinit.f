       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer i,j, blockno
      double precision xlow, ylow, w

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      double precision init_radius, init_location(2)
      common /initradius_comm/ init_radius, init_location

      double precision ravg, x0, y0, xp,yp, zp,xc,yc, r, r0,B

      blockno = 0

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(1,i,j) = w

c            xc = xlow + dx/2.d0
c            yc = ylow + dy/2.d0
c            call mapc2m_annulus(xc,yc,xp,yp,zp)
c
c            x0 = init_location(1)
c            y0 = init_location(2)
c
c            r0 = init_radius
c            B = 1 !! maximum value
c            r = B + (-B/r0**2)*((xp - x0)**2 + (yp-y0)**2)
c            if (r .gt. 0) then
c               q(1,i,j) = r
c            else
c               q(1,i,j) = 0
c            endif
         enddo
      enddo

      return
      end


