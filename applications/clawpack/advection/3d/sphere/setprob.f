      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer manifold
      common /com_manifold/ manifold

      double precision maxelev, revs_per_second
      common /com_sphere/ maxelev, revs_per_second

      integer example

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) manifold
      read(10,*) revs_per_second
      read(10,*) maxelev  !! Can we get this to the z mapping? 
      close(10)

      end


