      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_second
      common /spherecomm/ revs_per_second

      integer manifold
      common /com_manifold/ manifold

      double precision maxelev 
      common /com_sphere/ maxelev

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

      double precision function sphere_get_maxelev()
      implicit none

      double precision maxelev 
      common /com_sphere/ maxelev

      sphere_get_maxelev = maxelev

      return

      end



