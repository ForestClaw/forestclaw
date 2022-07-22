c     # ------------------------------------------------------------------
c     # MAPC2M_PILLOWDISK
c     # ------------------------------------------------------------------
c     #
c     # Maps a logically rectangular Cartesian grid in [-1,1]x[-1,1] to
c     # the unit disk.
c     #
c     # ------------------------------------------------------------------
      subroutine mapc2m_pillowdisk(blockno,xc1,yc1,xp,yp,zp)
      implicit none

      double precision xc1,yc1,xp,yp,zp
      double precision xc,yc,zc
      integer blockno

c     # Map to [-1,1]x[-1,1]
      call mapc2m_cart(blockno,xc1,yc1,xc,yc,zc)

c     # Get circle of radius sqrt(2.d0)
      call mapc2p_disk_circle(xc,yc,xp,yp)

c     # Scale to get unit circle.
      xp = xp/sqrt(2.d0)
      yp = yp/sqrt(2.d0)
      zp = 0

      return
      end


      subroutine mapc2p_disk_circle(x1,y1,xp,yp)
      implicit none

      double precision xc,yc,xp,yp, x1,y1
      double precision xi,eta,x,y, minxy,maxxy

      xc = x1
      yc = y1

      xi = min(abs(xc),abs(yc))
      eta = max(abs(xc),abs(yc))
      eta = max(eta,1.d-10)

      call map_disk_north_sector(xi,eta,x,y)

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

      subroutine map_disk_north_sector(xi,eta,x,y)
      implicit none

      double precision xi,eta,x,y

      x = xi
      y = sqrt(2 - xi**2) - sqrt(2 - eta**2) + eta

      end
