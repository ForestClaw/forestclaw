      subroutine setprob
      implicit none

      double precision gamma,gamma1
      common /param/  gamma,gamma1

      open(unit=7,file='setprob.data')

      read(7,*) gamma
      gamma1 = gamma - 1.d0

      return
      end
