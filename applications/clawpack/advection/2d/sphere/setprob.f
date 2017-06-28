      subroutine sphere_setprob(revs_per_second)
      implicit none

      double precision revs_per_second, rps

      double precision pi
      common /compi/ pi

      common /spherecomm/ rps

      pi = 4.d0*atan(1.d0)

      rps = revs_per_second

      end
