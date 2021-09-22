      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision uvel, vvel, revs_per_s
      common /comm_velocity/ uvel, vvel, revs_per_s

      integer example
      common /comex/ example

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) uvel
      read(10,*) vvel
      read(10,*) revs_per_s

      end
