      subroutine setprob_annulus()
      implicit none

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

      double precision revs_per_s, vertical_speed
      common /stream_comm/ revs_per_s, vertical_speed

c     # Must use edge velocities
      integer color_equation
      common /eqn_comm/ color_equation

c     # used only for edge velocities
      integer use_stream
      common /velocity_comm/ use_stream

      double precision beta
      common /annulus_comm/ beta

      integer refine_pattern
      common /refine_comm/ refine_pattern

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      iunit = 10
      fname = 'setprob.data'      
      open(iunit,file=fname)

      read(iunit,*) example
      read(iunit,*) mapping
      read(iunit,*) initchoice
      read(iunit,*) revs_per_s
      read(iunit,*) twist
      read(iunit,*) vertical_speed
      read(iunit,*) color_equation
      read(iunit,*) use_stream
      read(iunit,*) beta
      read(iunit,*) refine_pattern
      read(iunit,*) init_radius    !! radius
      close(iunit)

      open(iunit,file='twist.dat')
      if (mapping == 0) then
            write(iunit,*) 0.0
      else
            write(iunit,*) twist
      endif
      close(iunit)

      end
