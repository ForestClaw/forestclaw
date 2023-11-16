      subroutine fclaw_map_2d_c2m_identity(blockno,xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno

c     # Map in [0,1]x[0,1]
      xp = xc
      yp = yc
      zp = 0


      end
