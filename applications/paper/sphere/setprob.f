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

      integer refine_pattern
      common /refine_comm/ refine_pattern

      double precision omega(3)
      common /rotation_comm/ omega

      double precision x0(2), y0(2), z0(2)
      common /qinit_comm/ x0, y0, z0

      double precision b_init, c_init
      common /init_parms/ b_init, c_init

      double precision sharpness
      common /hsmooth_parms/ sharpness

      integer curvature_correction
      common /conservation_com/ curvature_correction


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
      read(10,*) b_init
      read(10,*) c_init
      read(10,*) sharpness
      read(10,*) curvature_correction;
      close(10)


c    # Center for two disks  
      if (example .eq. 0) then
          x0(1) = 0.866025403784439d0    !! cos(5*pi/6)
          x0(2) = 0.866025403784439d0    !! cos(7*pi/6)
          y0(1) = 0.5
          y0(2) = -0.5
      endif




      end
