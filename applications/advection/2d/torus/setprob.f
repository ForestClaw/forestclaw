      subroutine torus_setprob(example_in,alpha_in)
      implicit none

      double precision pi
      common /compi/ pi

      double precision alpha_in, alpha
      common /toruscomm/ alpha

      integer example_in, example
      common /excomm_example/ example

      pi = 4.d0*atan(1.d0)

      example = example_in
      alpha = alpha_in

      end
