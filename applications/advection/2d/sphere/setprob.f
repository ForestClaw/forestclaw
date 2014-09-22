      subroutine setprob()
      implicit none

      double precision kappa,tfinal

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # These are the values that work well with the hemisphere.  Other
c     # values don't work so well, since we don't have any inflow
c     # conditions set.

      kappa = 0
      tfinal = 5.0
      call set_wind_parms(kappa,tfinal);

      end
