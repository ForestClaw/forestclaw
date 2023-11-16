      subroutine fclaw_map_2d_c2m_pillowdisk5(blockno,xc1,yc1,
     &      xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc1, yc1
      double precision xc,yc,zc,xp,yp,zp,alpha

      call fclaw_map_2d_c2m_fivepatch(blockno,xc1,yc1,xc,yc,zc,alpha)

c     # Map from [-1,1]x[-1,1] back to [0,1]x[0,1]
      xc = (xc + 1)/2.d0
      yc = (yc + 1)/2.d0

      call fclaw_map_2d_c2m_pillowdisk(blockno,xc,yc,xp,yp,zp)

      end
