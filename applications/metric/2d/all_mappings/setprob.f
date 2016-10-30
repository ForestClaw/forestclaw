      subroutine metric_setprob(beta_in)
      implicit none

      double precision beta_in, beta
      common /comtorus/ beta

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

      beta = beta_in

      end
