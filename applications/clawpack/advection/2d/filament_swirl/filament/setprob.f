      subroutine filament_setprob
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer solver
      common /solver/ solver

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      solver = 1

      end
