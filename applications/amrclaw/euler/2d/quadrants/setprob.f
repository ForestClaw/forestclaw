      subroutine setprob
      implicit none

      double precision gamma,gamma1
      common /param/  gamma,gamma1

c     # This should be set up using options
      gamma = 1.4d0
      gamma1 = gamma-1.d0

      return
      end
