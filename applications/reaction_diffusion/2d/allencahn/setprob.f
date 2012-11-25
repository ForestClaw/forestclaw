      subroutine setprob()
      implicit none

      double precision pi
      common /compi/ pi

      double precision D
      common /comsrc/ D

      pi = 4.d0*atan(1.d0)

      open(20,file='setprob.data')

      read(20,*) D
      close(20)


      end
