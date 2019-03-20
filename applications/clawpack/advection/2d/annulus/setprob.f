      subroutine setprob_annulus(example_in,mapping_in, ic_in, 
     &     revs_per_s_in, ceqn_in, use_stream_in, beta_in,
     &     refine_pattern_in)
      implicit none


      integer example_in, mapping_in, ic_in, ceqn_in      
      integer use_stream_in, refine_pattern_in
      double precision beta_in, revs_per_s_in

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example

      double precision beta
      common /annulus_comm/ beta

      integer use_stream
      common /velocity_comm/ use_stream

      integer mapping
      common /mapping_comm/ mapping

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer initchoice
      common /initchoice_comm/ initchoice

      integer color_equation
      common /eqn_comm/ color_equation

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

      example = example_in
      mapping = mapping_in
      refine_pattern = refine_pattern_in

      beta = beta_in

      initchoice = ic_in

      revs_per_s = revs_per_s_in

      color_equation = ceqn_in
      use_stream = use_stream_in

      end
