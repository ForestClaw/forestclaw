      subroutine replicated_setprob(example_in)
      implicit none

      integer example_in, example

      double precision pi
      common /compi/ pi

      double precision u0_comm,v0_comm, revs_comm
      common /comm_velocity/ u0_comm,v0_comm,revs_comm
      common /comex/ example

      pi = 4.d0*atan(1.d0)

      example = example_in

      revs_comm = 0.5d0
      u0_comm = 1.d0
      v0_comm = 1.d0

      end
