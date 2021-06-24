      subroutine square_setprob()

      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision alpha
      common /fivepatch_comm/ alpha

      double precision x0, y0
      common /bilinear_comm/ x0, y0

      double precision velocity(2)
      common /velocity_comm/ velocity

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) mapping
      read(10,*) initchoice
      read(10,*) alpha
      read(10,*) x0
      read(10,*) y0
      read(10,*) velocity(1)
      read(10,*) velocity(2)
      close(10)


      end
