      subroutine setprob_torus(example)
      implicit none

      integer example

      double precision pi
      common /compi/ pi

      integer ex_comm
      common /comm_example/ ex_comm

      double precision u0_comm,v0_comm, revs_comm
      common /comm_velocity/ u0_comm,v0_comm,revs_comm

      double precision revs_per_sec_comm, r0_comm
      common /comm_torus/ revs_per_sec_comm, r0_comm

      pi = 4.d0*atan(1.d0)

      ex_comm = example

      if (example .eq. 6) then
         revs_comm = 0.5d0
         u0_comm = 1.d0
         v0_comm = 1.d0
      elseif (example .eq. 7) then
         revs_per_sec_comm = 0.5d0
         r0_comm = 0.4d0
      endif

      end
