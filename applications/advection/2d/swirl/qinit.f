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

       integer i, j
       double precision xi, yj, s

       do i = 1-mbc,mx+mbc
          xi = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy

             s = 79.d0/128.d0
             s = 0.5d0
c             s = 5.d0/8.d0 + 1.d0/128.d0
             if (xi .lt. s) then
                q(i,j,1) = 1.d0
             else
                q(i,j,1) = 0.d0
             endif
          enddo
       enddo

       return
       end
