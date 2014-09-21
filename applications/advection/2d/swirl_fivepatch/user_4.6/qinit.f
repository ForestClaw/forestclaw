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
       double precision xlow,ylow,w

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             xlow = xlower + (i-0.5d0)*dx
             do j = 1-mbc,my+mbc
                ylow = ylower + (j-1d0)*dy
                call cellave2(xlow,ylow,dx,dy,w)
                q(i,j,1) = w
             enddo
          enddo
       enddo

       return
       end

      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer*8 cont, get_context

      integer blockno, get_block
      double precision r

      logical fclaw2d_map_is_used

      cont = get_context()
      blockno = get_block()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)
      else
         xp = xc
         yp = yc
      endif

      fdisc = xp-0.5d0
      end
