      subroutine setprob_torus(example)
      implicit none

      double precision kappa,tfinal
      integer example, ex_comm

      double precision pi
      common /compi/ pi
      common /comm_example/ ex_comm

      pi = 4.d0*atan(1.d0)

      ex_comm = example

      end
