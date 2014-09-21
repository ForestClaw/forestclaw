c     =================================================
      double precision function psi(x,y)
c     =================================================

c     # stream function

      implicit none

      double precision x,y,pi
      common /compi/ pi

      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      return
      end
