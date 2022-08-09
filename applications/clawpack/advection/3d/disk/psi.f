      double precision function psi(x,y,z,t)
      implicit none

      double precision x,y,z,t, r

      double precision pi, pi2
      common /compi/ pi, pi2

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

c     # Rigid body rotation
c      psi = r**2

c     # Filament formation.  The filament is in a 2x2 box, so 
c     # we have to scale "dx" and "dy" by four.  
      psi = (4.d0/3.d0)*r**3

      return
      end