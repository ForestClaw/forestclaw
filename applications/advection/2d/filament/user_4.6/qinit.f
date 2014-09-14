c     =====================================================
       subroutine qinit_filament(meqn,mbc,mx,my,xlower,ylower,
     &      dx,dy,q,maux,aux,ismanifold)
c     =====================================================

c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise

       implicit none

       integer maxmx, maxmy, meqn, mbc, mx, my, maux
       double precision xlower, ylower, dx, dy
       double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
       double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
       logical ismanifold, ismanifold_com

       integer i, j, mq
       double precision xlow, ylow, w

       common /commanifold/ ismanifold_com

       ismanifold_com = ismanifold

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

      r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2)

      fdisc = r-0.25d0
c      fdisc = xp - 1
      end
