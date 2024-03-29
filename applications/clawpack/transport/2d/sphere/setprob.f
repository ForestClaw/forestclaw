      subroutine setprob()

      implicit none

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      integer initchoice
      common /initchoice_comm/ initchoice

      integer refine_pattern
      common /refine_comm/ refine_pattern

      double precision omega(3)
      common /rotation_comm/ omega

      double precision x0(2), y0(2), z0(2)
      common /qinit_comm/ x0, y0, z0

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      open(10,file='setprob.data')
      read(10,*) example
      read(10,*) mapping
      read(10,*) initchoice
      read(10,*) omega(1)
      read(10,*) omega(2)
      read(10,*) omega(3)
      read(10,*) refine_pattern
      close(10)

c    # Center for two disks  
      if (example .eq. 0) then
          x0(1) = 0.866025403784439d0    !! cos(5*pi/6)
          x0(2) = 0.866025403784439d0    !! cos(7*pi/6)
          y0(1) = 0.5
          y0(2) = -0.5
      endif




      end
