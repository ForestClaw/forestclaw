      subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      integer blockno, clawpack_get_block
      double precision x,y,z, xlow, ylow, w

      blockno = clawpack_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(i,j,1) = w
         enddo
      enddo

      return
      end


      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp
      logical iscart
      logical fclaw2d_map_is_used

      double precision pi
      common /compi/ pi

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)
         if (iscart()) then
            rp = sqrt(xp**2 + yp**2)
            fdisc = rp-0.25d0
         else
            th = atan2(yp,xp)
            tp = abs(th + pi/2.d0)
            fdisc = tp - pi/8.d0
         endif
      else
         xp = xc
         yp = yc
         rp = sqrt(xp**2 + yp**2)
         fdisc = rp-0.25d0
      endif

      end
