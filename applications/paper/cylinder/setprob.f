      subroutine cylinder_setprob()
      implicit none

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
      close(10)
      
c      h_cyl = 2*pi*r_cyl

      end
