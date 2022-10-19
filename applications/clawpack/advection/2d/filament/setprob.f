      subroutine setprob
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /com_ex/ example

      double precision alpha
      common /com_alpha/ alpha

      double precision center(2)
      common /com_bilinear/ center
      
      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) alpha
      read(10,*) center(1)
      read(10,*) center(2)
      close(10)

      end
