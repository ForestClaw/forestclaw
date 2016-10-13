      subroutine swirl_setprob(period_in)
      implicit none

      double precision pi, tperiod, pi2, period_in

      common /compi/ pi
      common /comvt/ tperiod,pi2

      tperiod = period_in

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      end
