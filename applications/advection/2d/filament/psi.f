c     =================================================
      double precision function psi(x,y)
c     =================================================

c     # stream function

      implicit none

      double precision x,y,pi,r
      common /compi/ pi

c      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)
      psi = (4.d0/3.d0)*r**3

c      psi = x - y

      return
      end
