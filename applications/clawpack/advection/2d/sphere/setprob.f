      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_second
      common /spherecomm/ revs_per_second

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) revs_per_second
      close(10)

      end
