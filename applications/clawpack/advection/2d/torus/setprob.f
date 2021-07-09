      subroutine setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      integer example
      common /example_comm/ example

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) alpha
      read(10,*) beta
      read(10,*) revs_per_s
      close(10)

      end
