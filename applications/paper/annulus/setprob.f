      subroutine setprob_annulus()
      implicit none

      integer iunit
      character(len=25) fname      

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      double precision revs_per_s, cart_speed, amplitude,freq
      common /velocity_comm/ revs_per_s, cart_speed, amplitude, freq

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer refine_pattern
      common /refine_comm/ refine_pattern

      double precision init_radius
      common /initradius_comm/ init_radius

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      iunit = 10
      fname = 'setprob.data'      
      open(iunit,file=fname)

      read(iunit,*) example
      read(iunit,*) revs_per_s
      read(iunit,*) cart_speed
      read(iunit,*) amplitude
      read(iunit,*) freq
      read(iunit,*) beta
      read(iunit,*) theta(1), theta(2)
      read(iunit,*) refine_pattern
      read(iunit,*) init_radius    !! radius
      close(iunit)

      cart_speed = 1.092505803290319d0


      open(10,file='mapping.dat')
      write(10,*) example
      write(10,*) amplitude
      write(10,*) init_radius
      write(10,*) beta
      write(10,*) theta(1)
      write(10,*) theta(2)
      write(10,*) freq
      write(10,*) cart_speed
      close(10)

      end
