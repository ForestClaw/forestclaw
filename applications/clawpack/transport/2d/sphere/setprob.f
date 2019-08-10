      subroutine sphere_setprob()

      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision x0(2), y0(2)
      common /qinit_comm/ x0, y0

      double precision kappa, period
      common /wind_comm/ kappa, period

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) mapping
      read(10,*) initchoice
      read(10,*) kappa
      read(10,*) period
      close(10)

c      x0(1) = 5*pi/6.d0
c      x0(2) = 7*pi/6.d0

       x0(1) = 0.866025403784439d0
       x0(2) = 0.866025403784439d0

       y0(1) = 0.5
       y0(2) = -0.5



      end
