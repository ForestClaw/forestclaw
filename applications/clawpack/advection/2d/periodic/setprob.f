      subroutine periodic_setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi,pi2

      double precision uvel, vvel
      common /comvelocity/ uvel, vvel

      open(10,file='setprob.data')
      read(10,*) uvel
      read(10,*) vvel
      close(10)

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      end
