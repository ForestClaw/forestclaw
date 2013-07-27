c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================

c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise

       implicit none

       integer maxmx, maxmy, meqn, mbc, mx, my, maux
       double precision xlower, ylower, dx, dy
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       integer i, j,mq
       double precision xlow, ylow, w

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             do j = 1-mbc,my+mbc

                xlow = xlower + (i-1)*dx
                ylow = ylower + (j-1)*dy
                call cellave2(xlow,ylow,dx,dy,w)
                q(i,j,mq) = w

             enddo
          enddo
       enddo

       return
       end

      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc
      double precision xs,ys, rs

      xs = 0.25d0
      ys = 0.25d0

      rs = sqrt((xc - xs)**2 + (yc - ys)**2)
      fdisc = rs - 0.2d0
      end
