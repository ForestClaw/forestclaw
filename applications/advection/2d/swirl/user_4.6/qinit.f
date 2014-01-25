c     =====================================================
       subroutine qinit(maxmx,maxmy, meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================

c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise

       implicit none

       integer meqn, mbc, mx, my, maux, maxmx, maxmy
       double precision xlower, ylower, dx, dy
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       integer i, j, mq
       double precision xi, yj, s

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             xi = xlower + (i-0.5d0)*dx
             do j = 1-mbc,my+mbc
                yj = ylower + (j-0.5d0)*dy

                s = 0.5d0
                if (xi .lt. s) then
                   q(i,j,mq) = 0.d0
                else
                   q(i,j,mq) = 1.d0
                endif
             enddo
          enddo
       enddo

       return
       end
