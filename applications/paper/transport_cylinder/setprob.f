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

      double precision r_cyl, h_cyl
      common /cylinder_comm/ r_cyl, h_cyl

      integer exact_metric
      common /metric_comm/ exact_metric

      double precision xc0, yc0, r0
      common /cylinder_init_comm/ xc0, yc0, r0

      double precision revs_per_s, v_speed
      common /stream_comm/ revs_per_s, v_speed

      integer fluctuation_sync
      common /fluctuation_com/ fluctuation_sync

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
      read(10,*) r_cyl
      read(10,*) h_cyl
      read(10,*) xc0
      read(10,*) yc0
      read(10,*) r0
      read(10,*) revs_per_s
      read(10,*) v_speed
      read(10,*) fluctuation_sync
      read(10,*) exact_metric


      read(iunit,*) grid_mx
      read(iunit,*) mi
      read(iunit,*) mj
      read(iunit,*) maxlevel
      read(iunit,*) rfactor
      read(iunit,*) qad_mode
      close(iunit)

c     # Write out data to be read back into Matlab.  We could just 
c     # read in variables from setprob.data, but the Python scripts
c     # add extra comments that make it a pain to read. 

      open(10,file='mapping.dat')
      write(10,*) example
      write(10,*) initchoice
      write(10,*) refine_pattern
      write(10,*) r_cyl
      write(10,*) h_cyl
      write(10,*) xc0
      write(10,*) yc0
      write(10,*) r0
      write(10,*) revs_per_s
      write(10,*) v_speed
      close(10)

      end
