      subroutine setprob()

      implicit none

      integer iunit
      character(len=25) fname      

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

      double precision revs_per_s, cart_speed
      common /stream_comm/ revs_per_s, cart_speed

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      integer maxlevel, rfactor, grid_mx, mi, mj
      common /amr_comm/ maxlevel, rfactor, grid_mx, mi, mj

      integer qad_mode
      common /qad_comm/ qad_mode      


      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      iunit = 10
      fname = 'setprob.data'      
      call opendatafile(iunit, fname)
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


      read(iunit,*) grid_mx
      read(iunit,*) mi
      read(iunit,*) mj
      read(iunit,*) maxlevel
      read(iunit,*) rfactor
      read(iunit,*) qad_mode
      close(iunit)

      open(10,file='mapping.dat')
      write(10,*) example
      write(10,*) initchoice
      write(10,*) refine_pattern
      write(10,*) alpha
      write(10,*) beta
      write(10,*) init_radius
      write(10,*) revs_per_s
      write(10,*) cart_speed
      write(10,*) theta_range(1)
      write(10,*) theta_range(2)
      write(10,*) phi_range(1)
      write(10,*) phi_range(2)
      close(10)

      end
