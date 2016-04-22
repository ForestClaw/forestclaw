c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================

c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise

       implicit none

       integer meqn, mbc, mx, my, maux
       double precision xlower, ylower, dx, dy
       double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

       integer i, j, mq
       double precision xi, yj, s
       double precision xc,yc,zc,xp,yp,zp

       integer*8 cont, get_context
       integer blockno, fc2d_clawpack5_get_block

       blockno = fc2d_clawpack5_get_block()
       cont = get_context()

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             xi = xlower + (i-0.5d0)*dx
             do j = 1-mbc,my+mbc
                call fclaw2d_map_c2m(cont,
     &                blockno,xi,yj,xp,yp,zp)
                yj = ylower + (j-0.5d0)*dy

                s = 0.5d0
                if (xp .lt. s) then
                   q(mq,i,j) = 0.d0
                else
                   q(mq,i,j) = 1.d0
                endif
             enddo
          enddo
       enddo

       return
       end
