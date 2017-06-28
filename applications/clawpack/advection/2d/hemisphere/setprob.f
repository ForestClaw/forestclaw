      subroutine setprob()
      implicit none

      double precision pi
      common /compi/ pi

      double precision revs_per_second
      common /spherecomm/ revs_per_second

      pi = 4.d0*atan(1.d0)

      revs_per_second = 0.5d0

      end
