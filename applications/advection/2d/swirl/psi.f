c     =================================================
      double precision function psi(x,y)
c     =================================================

c     # stream function

      implicit none

      double precision x,y,pi
      common /compi/ pi

      psi = ((dsin(pi*x))**2 * (dsin(pi*y))**2) / pi

      return
      end
