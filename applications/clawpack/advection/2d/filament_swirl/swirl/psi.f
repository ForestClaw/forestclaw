      double precision function swirl_psi(xp,yp)
      implicit none

      double precision xp,yp

      double precision pi, pi2
      common /compi/ pi, pi2

      swirl_psi = ((sin(pi*xp))**2 * (sin(pi*yp))**2) / pi

      return
      end
