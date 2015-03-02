      subroutine quadrants_setprob(gamma)
      implicit none

      double precision gamma,gamma_com,gamma1
      common /param/  gamma_com,gamma1

c     # This should be set up uesing options
      gamma_com = gamma
      gamma1 = gamma-1.d0

      return
      end
