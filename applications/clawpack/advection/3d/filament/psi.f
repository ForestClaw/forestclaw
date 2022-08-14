      double precision function psi(x,y,z)
      implicit none

      double precision x,y,z,pi,r
      common /compi/ pi

c      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

c     # Filament formation (negative for clockwise rotation)
      psi = (4.d0/3.d0)*r**3

      return
      end
