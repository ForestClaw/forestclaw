      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
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

      integer example
      common /comm_example/ example

      blockno = fc2d_clawpack46_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            if (example .eq. 6) then
c              # use this to get a smooth solution for computing the error
               q(i,j,1) = q0(blockno,xc,yc)
            else
               call cellave2(blockno,xlow,ylow,dx,dy,w)
               q(i,j,1) = w
            endif
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
      double precision xloc(0:4),yloc(0:4),zloc(0:4)
      logical iscart
      logical fclaw2d_map_is_used
      integer i, mi, mj
      double precision xc1, yc1, zc1

      data xloc /0, 1, 1, 0, 0.5d0/
      data yloc /0, 0, 1, 1, 0.5d0/

      double precision pi
      common /compi/ pi

      integer example
      common /comm_example/ example

      cont = get_context()

      r0 = 0.4d0

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)

         if (iscart()) then
            do i = 0,4
               rp = sqrt((xp-xloc(i))**2 + (yp-yloc(i))**2)
               fdisc = rp-0.3d0
               if (fdisc .lt. 0) then
                  return
               endif
            enddo
         else
            do i = 0,4
               th = 2*pi*i/5
               xloc(i) = cos(th)
               yloc(i) = sin(th)
               zloc(i) = r0
            enddo

c           # Torus or annulus
            th = atan2(yp,xp)
            tp = abs(th + pi/2.d0)
            fdisc = tp - pi/8.d0

c           # Sphere
            do i = 0,4
               r2 = (xp - xloc(i))**2 + (yp-yloc(i))**2 +
     &               (zp-zloc(i))**2 - r0**2
               if (r2 < 0) then
                  fdisc = r2
                  return
               endif
            enddo
            fdisc = 1
         endif
      else
         if (example == 5) then
c           # Map each brick to a [0,1]x[0,1] domain and duplicate
c           # initial conditions.
            xp = xc
            yp = yc
            do i = 0,4
               rp = sqrt((xp-xloc(i))**2 + (yp-yloc(i))**2)
               fdisc = rp-0.3d0
               if (fdisc .lt. 0) then
                  return
               endif
            enddo
         else
c           # No mapping.
            rp = (xc-0.5d0)**2 + (yc-0.5d0)**2
            fdisc = rp-(0.1d0)**2
         endif
      endif

      end

      double precision function  q0(blockno,xc1,yc1)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      double precision xc1, yc1
      integer blockno
      integer*8 cont, get_context
      double precision r,r0
      logical fclaw2d_map_is_used
      double precision Hsmooth,h1,h2
      integer i

      double precision xloc(0:4),yloc(0:4),zloc(0:4)

      double precision pi
      common /compi/ pi

      integer example
      common /comm_example/ example

      data xloc /0, 1, 1, 0, 0.5d0/
      data yloc /0, 0, 1, 1, 0.5d0/

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         xc = xc1
         yc = yc1
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)

         r0 = 0.4d0
         r = sqrt((xp - 1.0)**2 + yp**2 + (zp-r0)**2)
         q0 = Hsmooth(r + r0) - Hsmooth(r - r0)
      else
         if (example .eq. 6) then
            if (xc1 .lt. 0) then
               xp = 1 - (-xc1 - int(-xc1))
            else
               xp = xc1 - int(xc1)
            endif
            if (yc1 .lt. 0) then
               yp = 1-(-yc1 - int(-yc1))
            else
               yp = yc1 - int(yc1)
            endif
            r0 = 0.2d0
            rp = sqrt((xp-xloc(4))**2 + (yp-yloc(4))**2)
            h1 = Hsmooth(rp+r0)
            h2 = Hsmooth(rp-r0)
            q0 = h1 - h2
         else
            xp = xc
            yp = yc
            rp = sqrt(xp**2 + yp**2)
            q0 = rp-0.25d0
         endif
      endif

      end

      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end
