      subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      integer blockno, fc2d_clawpack46_get_block
      double precision x,y,z, xlow, ylow, w,xc,yc
      double precision q0

      blockno = fc2d_clawpack46_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(i,j,1) = w

c            xc = xlow + dx/2.d0
c            yc = ylow + dy/2.d0
c            q(i,j,1) = q0(blockno,xc,yc)
         enddo
      enddo

      return
      end


      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp, r2, r0
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
c           # Torus or annulus
            th = atan2(yp,xp)
            tp = abs(th + pi/2.d0)
            fdisc = tp - pi/8.d0

c           # Sphere
            r0 = 0.4d0
    2       r2 = (xp - 1.0)**2 + yp**2 + (zp-0.4d0)**2
            fdisc = r2 - r0**2
         endif
      else
         xp = xc
         yp = yc
         rp = sqrt(xp**2 + yp**2)
         fdisc = rp-0.25d0
      endif

      end

      double precision function  q0(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno
      integer*8 cont, get_context
      double precision r,r0
      logical fclaw2d_map_is_used
      double precision Hsmooth

      double precision pi
      common /compi/ pi

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)

         r0 = 0.4d0
         r = sqrt((xp - 1.0)**2 + yp**2 + (zp-r0)**2)
         q0 = Hsmooth(r + r0) - Hsmooth(r - r0)
      else
         xp = xc
         yp = yc
         rp = sqrt(xp**2 + yp**2)
         q0 = rp-0.25d0
      endif

      end

      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.05d0) + 1)/2.0

      end
