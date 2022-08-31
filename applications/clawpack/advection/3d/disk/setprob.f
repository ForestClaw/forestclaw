      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /com_ex/ example

      double precision alpha
      common /com_alpha/ alpha

      integer manifold
      common /com_manifold/ manifold
      
      pi = 4.d0*atan(1.d0)
      pi2 = 2.0*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) manifold  !!  Needed for 3d Riemann solvers for advection
      read(10,*) alpha
      close(10)


      end
