      subroutine setprob()
      implicit none

      double precision kappa,tfinal

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

      kappa = 0
      tfinal = 5.0
      call set_wind_parms(kappa,tfinal);

      end
