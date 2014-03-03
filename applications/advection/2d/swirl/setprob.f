      subroutine setprob
      implicit none

      double precision pi, tperiod, pi2

      common /compi/ pi
      common /comvt/ tperiod,pi2

c      open(unit=7,file='setprob.data')
c      read(7,*) tperiod
c      close(7)
      tperiod = 4.d0

      call set_maptype_cart()

      pi = 4.d0*atan(1.d0)
      pi2 = 2.d0*pi


      end
