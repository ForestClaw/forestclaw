      subroutine clawpack46_qinit(maxmx,maxmy, maxmz, meqn,mbc,
     &      mx,my,mz, xlower,ylower,zlower, dx,dy,dz, q,maux,aux)
       implicit none

       integer meqn, mbc, mx, my, mz, maux, maxmx, maxmy, maxmz
       double precision xlower, ylower, zlower, dx, dy, dz
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       integer i, j, mq
       double precision xc,yc

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             xc = xlower + (i-0.5d0)*dx
             do j = 1-mbc,my+mbc
                yc = ylower + (j-0.5d0)*dy
                if (xc .lt. 0.5d0) then
                   q(i,j,mq) = 1.d0
                else
                   q(i,j,mq) = 0
                endif
             enddo
          enddo
       enddo

       return
       end
