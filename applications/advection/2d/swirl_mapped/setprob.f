      subroutine swirl_setprob(tperiod_in)
      implicit none

      double precision pi, tperiod_in, tperiod, pi2

      common /compi/ pi
      common /comvt/ tperiod,pi2

      tperiod = tperiod_in
      pi = 4.d0*atan(1.d0)
      pi2 = 2.d0*pi


      end
