      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      pi = 4.d0*atan(1.d0)
      pi2 = 2.0*pi

      end
