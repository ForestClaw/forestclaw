      subroutine setprob()
      implicit none

      double precision grav
      common /cparam/grav    !# gravitational parameter

      double precision x0, y0, r0
      common/cdisc/ x0,y0,r0

      double precision hin, hout
      common /comic/ hin,hout

      integer example
      common /comex/ example

      open(10,file='setprob.data')

      read(10,*) example
      read(10,*) grav
      read(10,*) x0
      read(10,*) y0
      read(10,*) r0
      read(10,*) hin
      read(10,*) hout
      
      close(10)

      return
      end
