      subroutine setprob
      implicit none

      double precision pi, tperiod, pi2

      common /compi/ pi
      common /comvt/ tperiod,pi2

      tperiod = 4.d0

      pi = 4.d0*atan(1.d0)
      pi2 = 2.d0*pi


      end
