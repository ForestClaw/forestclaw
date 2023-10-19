      subroutine setprob
      implicit none

      double precision pi
      common /compi/ pi

      integer manifold
      common /com_manifold/ manifold

      integer example

      pi = 4.d0*atan(1.d0)

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) manifold
      !! Other values in setprob.data aren't needed by fortran routines
      close(10)


      end
