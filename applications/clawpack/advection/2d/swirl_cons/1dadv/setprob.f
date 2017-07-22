      subroutine setprob
      implicit none

      double precision pi
      integer example

      common /compi/ pi
      common /comex/ example

      open(10,file='setprob.data')
      read(10,*) example
      close(10)

      pi = 4.d0*atan(1.d0)

      return
      end
