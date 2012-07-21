      subroutine setprob
      implicit none

      double precision pi, tperiod, pi2

      common /compi/ pi
      common /comvt/ tperiod,pi2

      open(unit=7,file='setprob.data',status='old',form='formatted')
      read(7,*) tperiod
      close(7)

      pi = 4.d0 * datan(1.d0)
      pi2 = 2.d0*pi


      end
