c     # ------------------------------------------------------------------
c     # MAPC2M_SPHERE
c     # ------------------------------------------------------------------
c     #
c     # Maps a logically rectangular Cartesian grid in [-3,1]x[-1,1] to
c     # the unit sphere where [-1,1]x[-1,1] is mapped to the upper
c     # hemisphere
c     #
c     # ------------------------------------------------------------------
      subroutine mapc2m_sphere(xc,yc,xp,yp,zp)
      implicit none

      double precision x1,y1,xp,yp,zp
      double precision xc,yc, x2
      double precision tol, rp2, rp

      integer blockno, get_block


c     # Map ghost cells to interior of the domain
      tol = 0
      if (yc + tol .gt. 1.d0) then
         y1 = -yc + 2.d0
      elseif (yc - tol .lt. -1.d0) then
         y1 = -yc - 2.d0
      else
         y1 = yc
      endif

      if (xc + tol .lt. -3.d0) then
         if (abs(yc)+tol .gt. 1.d0) then
            x1 = -xc - 6
         else
            x1 = xc + 4
         endif
      elseif (xc - tol .gt. 1.d0) then
         if (abs(yc) + tol .gt. 1.d0) then
            x1 = -xc + 2
         else
             x1 = xc - 4.d0
         endif
      else
         if (abs(yc) + tol .gt. 1) then
            x1 = -xc - 2
         else
            x1 = xc
         endif
      endif

c     # Map lower hemisphere to upper hemisphere
      if (x1+tol .lt. -1.d0) then
         x2 = -(x1 + 2.d0)
      else
         x2 = x1
      endif

c     # Get circle of radius sqrt(2.d0)
      call mapc2p_circle_sp(x2,y1,xp,yp)

      rp2 = xp**2 + yp**2
      if (abs(rp2 - 2.d0) .lt. 1d-10) then
         zp = 0.d0
      else
         zp = sqrt(2.d0 - rp2)
      endif

c      zp = sqrt(2.d0-min(rp2,2.d0))

      if (x1 .lt. -1.d0) then
         zp = -zp
      endif

c     # This maps everything to the unit sphere - no vertical
c     # distance
      xp = xp/sqrt(2.d0)
      yp = yp/sqrt(2.d0)
      zp = zp/sqrt(2.d0)

      blockno = get_block()
      if (blockno .eq. 1) then
         zp = -zp
      endif

      return
      end
c --------------------------------- MAPC2M -------------------------------


      subroutine mapc2p_circle_sp(x1,y1,xp,yp)
      implicit none

      double precision xc,yc,xp,yp, x1,y1
      double precision xi,eta,x,y, minxy,maxxy
      double precision xit, etat, dd

      double precision pi
      common /compi/ pi

      xc = x1
      yc = y1

      if (x1 .lt. -1.d0) then
         xc = -(xc + 2)
      endif

      xi = min(abs(xc),abs(yc))
      eta = max(abs(xc),abs(yc))
      eta = max(eta,1.d-10)

      dd = sin(pi*eta/2.d0)
c      dd = eta*(2-eta)

      xit = (xi/eta)*dd
      etat = dd

      call map_north_sector_sp(xit,etat,x,y)

      minxy = min(abs(x),abs(y))
      maxxy = max(abs(x),abs(y))

      if (abs(xc) .le. abs(yc))  then
c        # In NS sectors
         xp = sign(1.d0,xc)*minxy
         yp = sign(1.d0,yc)*maxxy
      else
c        # In EW sectors
         xp = sign(1.d0,xc)*maxxy
         yp = sign(1.d0,yc)*minxy
      endif


      end

      subroutine map_north_sector_sp(xi,eta,x,y)
      implicit none

      double precision xi,eta,x,y

      x = xi
      y = sqrt(2 - xi**2) - sqrt(2 - eta**2) + eta

      end
