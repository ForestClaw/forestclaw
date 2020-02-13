      subroutine torus_setprob()
      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer color_equation
      common /eqn_comm/ color_equation

      integer use_stream
      common /velocity_comm/ use_stream

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) mapping
      read(10,*) initchoice
      read(10,*) alpha
      read(10,*) beta
      read(10,*) revs_per_s
      read(10,*) color_equation
      read(10,*) use_stream
      close(10)

      end
