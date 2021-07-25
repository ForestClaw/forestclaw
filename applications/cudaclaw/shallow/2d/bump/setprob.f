      subroutine bump_setprob()
      implicit none

      integer example
      common /comm_example/ example

      double precision grav
      common /cparam/grav    !# gravitational parameter

      open(10,file='setprob.data')
      read(10,*) example      
      read(10,*) grav
      close(10)


      return
      end
