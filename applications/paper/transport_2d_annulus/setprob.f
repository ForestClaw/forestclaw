      subroutine setprob()
      implicit none

      integer iunit
      character(len=25) fname      

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      double precision init_radius, init_location(2)
      common /initradius_comm/ init_radius, init_location

      double precision revs_per_s, vcart(2), amplitude, freq
      common /stream_comm/ revs_per_s, vcart, amplitude, freq

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer maxlevel, rfactor, grid_mx, mi, mj
      common /amr_comm/ maxlevel, rfactor, grid_mx, mi, mj

      integer refine_pattern, rotate_position
      common /refine_comm/ refine_pattern, rotate_position

      integer qad_mode
      common /qad_comm/ qad_mode      

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      iunit = 10
      fname = 'setprob.data'      
      call opendatafile(iunit, fname)

      read(iunit,*) example
      read(iunit,*) revs_per_s
      read(iunit,*) vcart(1)
      read(iunit,*) vcart(2)
      read(iunit,*) amplitude
      read(iunit,*) freq
      read(iunit,*) init_radius    !! radius
      read(iunit,*) init_location(1)
      read(iunit,*) init_location(2)
      read(iunit,*) beta
      read(iunit,*) theta(1)
      read(iunit,*) theta(2)

      read(iunit,*) grid_mx
      read(iunit,*) mi
      read(iunit,*) mj
      read(iunit,*) maxlevel
      read(iunit,*) rfactor
      read(iunit,*) refine_pattern
      read(iunit,*) rotate_position
      read(iunit,*) qad_mode
      close(iunit)

      open(10,file='mapping.dat')
      if (example .le. 1) then
            write(10,*) 0
      else
            write(10,*) amplitude
      endif
      write(10,*) init_radius
      write(10,*) init_location(1)
      write(10,*) init_location(2)
      write(10,*) beta
      write(10,*) theta(1)
      write(10,*) theta(2)
      write(10,*) vcart(1)
      write(10,*) vcart(2)
      write(10,*) freq
      close(10)

      end
