      subroutine cylinder_setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer initchoice
      common /initchoice_comm/ initchoice

      integer refine_pattern
      common /refine_comm/ refine_pattern

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /cylinder_comm/ alpha, beta, theta_range, phi_range

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
      read(10,*) init_radius
      read(10,*) revs_per_s
      read(10,*) cart_speed
      read(10,*) theta_range(1)
      read(10,*) theta_range(2)
      read(10,*) phi_range(1)
      read(10,*) phi_range(2)
      close(10)

      end
