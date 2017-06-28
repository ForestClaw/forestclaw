      subroutine setprob_annulus(beta_in)
      implicit none

      double precision pi
      common /compi/ pi

      double precision beta_in, beta
      common /commannulus/ beta

      pi = 4.d0*atan(1.d0)

      beta = beta_in

      end
