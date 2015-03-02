      subroutine swirl_setprob(tperiod)
      implicit none

      double precision pi, tperiod, pi2
      double precision tperiod_com

      common /compi/ pi
      common /comvt/ tperiod_com,pi2

      tperiod_com = tperiod

      pi = 4.d0*atan(1.d0)
      pi2 = 2.d0*pi


      end
