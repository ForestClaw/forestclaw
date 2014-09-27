      subroutine qinit_manifold(meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux,this_block_idx
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      double precision x,y,z, xlow, ylow, w

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(xlow,ylow,dx,dy,w)
            q(i,j,1) = w
         enddo
      enddo

      return
      end


      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer*8 cont, get_context
      integer blockno, get_block
      double precision th, tp
      logical iscart

      double precision pi
      common /compi/ pi

      cont = get_context()
      blockno = get_block()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      if (iscart()) then
         rp = sqrt(xp**2 + yp**2)
         fdisc = rp-0.25d0
      else
         th = atan2(yp,xp)
         tp = abs(th - pi/2.d0)
         fdisc = tp - pi/8.d0
      endif

      end
