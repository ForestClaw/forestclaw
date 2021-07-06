      subroutine cylinder_setprob()
      implicit none

      double precision arc

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer initchoice
      common /initchoice_comm/ initchoice

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer mapping 
      common /mapping_comm/ mapping

      double precision r_cyl, h_cyl
      common /cylinder_comm/ r_cyl, h_cyl

      double precision r_latlong, phi0, phi1
      common /latlong_comm/ r_latlong, phi0, phi1

      integer exact_metric
      common /metric_comm/ exact_metric

      double precision xc0, yc0, r0
      common /cylinder_init_comm/ xc0, yc0, r0

      double precision revs_per_s, v_speed
      common /stream_comm/ revs_per_s, v_speed

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

c     !! These are written out in cylinder_user.cpp
      open(10,file='setprob.data')
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
      read(10,*) exact_metric
      read(10,*) mapping
      close(10)
      
      r_latlong = sqrt((h_cyl/2)**2 + r_cyl**2)
      arc = asin(h_cyl/(2*r_latlong))
c     # Want phi = phi0 + (phi1-phi0)*y for y in [0,1]
c     # and phi in [-pi/2, pi/2]
      phi0 = -0.5*arc/(pi/2);
      phi1 = -phi0;


c      h_cyl = 2*pi*r_cyl

      end
