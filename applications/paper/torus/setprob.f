      subroutine torus_setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer initchoice
      common /initchoice_comm/ initchoice

      integer refine_pattern
      common /refine_comm/ refine_pattern

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision revs_per_s, cart_speed
      common /stream_comm/ revs_per_s, cart_speed

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) initchoice
      read(10,*) refine_pattern
      read(10,*) alpha
      read(10,*) beta
      read(10,*) revs_per_s
      read(10,*) cart_speed
      close(10)

      end
