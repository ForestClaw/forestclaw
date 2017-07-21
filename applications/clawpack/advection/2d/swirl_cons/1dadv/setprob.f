      subroutine setprob
      implicit double precision (a-h,o-z)

      double precision pi

      common /compi/ pi

      pi = 4.d0*atan(1.d0)

      return
      end
