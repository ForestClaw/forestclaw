      subroutine setprob()
      implicit none

      double precision gamma, gamma1
      common /cparam/  gamma, gamma1

      open(10,file='setprob.data')
      read(10,*) gamma
      close(10)

      gamma1 = gamma-1.d0

      return
      end
