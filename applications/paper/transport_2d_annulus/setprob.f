      subroutine setprob()
      implicit none

c     example_in,
c     mapping_in, 
c     ic_in, 
c     revs_per_s_in, 
c     ceqn_in, 
c     use_stream_in, 
c     beta_in


      integer iunit
      character(len=25) fname      

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      double precision twist
      common /twist_comm/ twist

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision revs_per_s, cart_speed
      common /stream_comm/ revs_per_s, cart_speed

c     # Must use edge velocities
      integer color_equation
      common /eqn_comm/ color_equation

c     # used only for edge velocities
      integer use_stream
      common /velocity_comm/ use_stream

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer maxlevel, rfactor, grid_mx, mi, mj
      common /amr_comm/ maxlevel, rfactor, grid_mx, mi, mj

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer qad_mode
      common /qad_comm/ qad_mode      

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      iunit = 10
      fname = 'setprob.data'      
      call opendatafile(iunit, fname)

      read(iunit,*) example
      read(iunit,*) mapping
      read(iunit,*) initchoice
      read(iunit,*) revs_per_s
      read(iunit,*) twist
      read(iunit,*) cart_speed
      read(iunit,*) init_radius    !! radius
      read(iunit,*) color_equation
      read(iunit,*) use_stream
      read(iunit,*) beta
      read(iunit,*) theta(1)
      read(iunit,*) theta(2)


      read(iunit,*) grid_mx
      read(iunit,*) mi
      read(iunit,*) mj
      read(iunit,*) maxlevel
      read(iunit,*) rfactor
      read(iunit,*) refine_pattern
      read(iunit,*) qad_mode
      close(iunit)

      end
