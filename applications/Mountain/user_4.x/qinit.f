       subroutine qinit(maxmx,maxmy, meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
       implicit none

       integer meqn, mbc, mx, my, maux, maxmx, maxmy
       double precision xlower, ylower, dx, dy
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       integer i, j, mq

       integer*8 cont, get_context
       integer blockno, fc2d_clawpack46_get_block
       double precision xc,yc,xp,yp,zp

       cont = get_context()
       blockno = fc2d_clawpack46_get_block()


       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             do j = 1-mbc,my+mbc
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
     &
                if (xp .lt. 1000) then
                   q(i,j,mq) = 1
                else
                   q(i,j,mq) = 0
                endif
             enddo
          enddo
       enddo

       return
       end
