      subroutine quadrants_setprob(gamma_in)
      implicit none

      double precision gamma_in,gamma,gamma1
      common /cparam/  gamma,gamma1

      gamma = gamma_in
      gamma1 = gamma-1.d0

      return
      end
