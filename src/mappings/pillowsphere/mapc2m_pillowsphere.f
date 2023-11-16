c     # ------------------------------------------------------------------
c     # MAPC2M_PILLOWSPHERE
c     # ------------------------------------------------------------------
c     #
c     # Maps a logically rectangular Cartesian grid in [-1,1]x[-1,1] to
c     # either the upper hemisphere (blockno == 0) or the lower hemisphere
c     # (blockno == 1).
c     #
c     # ------------------------------------------------------------------
      subroutine fclaw_map_2d_c2m_pillowsphere(
     &      blockno, xc1,yc1,xp,yp,zp)
      implicit none

      external fclaw_map_2d_c2m_cart
      double precision xc1,yc1, xp, yp, zp

      double precision x1, y1, d, rp2, xc, yc, zc
      integer blockno
      logical l1, l2
      logical ispillowsphere

c     # Map to [-1,1]x[-1,1]
      call fclaw_map_2d_c2m_cart(blockno,xc1,yc1,xc,yc,zc)

c     # Map xc and yc from ghost cells to interior
      d = max(xc - 1.0,0.d0) + max(-1 - xc,0.d0)
      x1 = ((1-d)/(1+d))*xc

      d = max(yc - 1,0.d0) + max(-1 - yc,0.d0)
      y1 = ((1-d)/(1+d))*yc

c     # Get circle of radius sqrt(2.d0).  Cluster radial
c     # direction towards boundary
      call mapc2p_circle_sp(x1,y1,xp,yp)

c     # Set z value
      rp2 = xp**2 + yp**2
      if (abs(rp2 - 2.d0) .lt. 1d-10) then
         zp = 0.d0
      else
         zp = sqrt(2.d0 - rp2)
      endif

c     # Values that are outside of [-1,1]x[-1,1] (in ghost cell regions)
c     # to the lower hemisphere
c     # Juqueen complained about the 'xor', because it is not in the
c     # Fortran standard.  Use .neqv. instead.
      l1 = abs(xc) .gt. 1
      l2 = abs(yc) .gt. 1
      if (l1 .neqv. l2) then
         zp = -zp
      endif

c     # This maps everything to the unit sphere
      xp = xp/sqrt(2.d0)
      yp = yp/sqrt(2.d0)
      zp = zp/sqrt(2.d0)

c     # Set lower hemisphere
      if (ispillowsphere()) then
         if (blockno .eq. 1) then
            zp = -zp
         endif
      endif

      return
      end
c --------------------------------- MAPC2M -------------------------------


c     # Map single grid to the disk.  Since this disk will be used for the
c     #
      subroutine mapc2p_circle_sp(x1,y1,xp,yp)
      implicit none

      double precision xc,yc,xp,yp, x1,y1
      double precision xi,eta,x,y, minxy,maxxy
      double precision xit, etat, dd

      double precision pi, pi2
      common /compi/ pi, pi2

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
