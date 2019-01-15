      subroutine torus_setprob(example_in,mapping_in, ic_in, 
     &     alpha_in,rps_in, ceqn_in)
      implicit none

      double precision alpha_in, rps_in
      integer example_in, mapping_in, ic_in, ceqn_in

      double precision pi
      common /compi/ pi

      double precision alpha
      common /torus_comm/ alpha

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example

      integer mapping
      common /mapping_comm/ mapping

      integer initchoice
      common /initchoice_comm/ initchoice

      integer color_equation
      common /eqn_comm/ color_equation

      pi = 4.d0*atan(1.d0)

      example = example_in
      mapping = mapping_in
      initchoice = ic_in

      alpha = alpha_in
      
      revs_per_s = rps_in

      color_equation = ceqn_in

      end
