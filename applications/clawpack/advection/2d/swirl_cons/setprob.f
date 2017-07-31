      subroutine swirl_setprob(ex_in)
      implicit none

      double precision pi, pi2
      integer example, ex_in

      common /compi/ pi
      common /comex/ example

      example = ex_in

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      end
