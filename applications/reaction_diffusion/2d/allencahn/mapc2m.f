c     # ------------------------------------------------------------------
c     # MAPC2M
c     # ------------------------------------------------------------------
c     #
c     # Mapping routine is defined here.  Main routine that is needed is
c     #
c     #        subroutine mapc2m(xc,yc,xp,yp,zp)
c     #
c     # ------------------------------------------------------------------
      subroutine mapc2m(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      xp = -1 + 2.d0*xc
      yp = -1 + 2.d0*yc
      zp = 0

      return
      end
