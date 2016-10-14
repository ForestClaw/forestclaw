      double precision function psi(xc,yc)
      implicit none

      double precision xc, yc
      double precision revs_per_s

      revs_per_s = 0.5d0

      psi = revs_per_s*(xc - yc)

      end
