      subroutine swirl_setprob(period_in, ex_in)
      implicit none

      double precision pi, tperiod, pi2, period_in
      integer example, ex_in

      common /compi/ pi
      common /comvt/ tperiod,pi2
      common /comex/ example

      example = ex_in

      tperiod = period_in

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      end
