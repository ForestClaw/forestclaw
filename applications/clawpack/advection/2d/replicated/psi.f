      double precision function psi(xc,yc)
      implicit none

      double precision xc, yc
      double precision revs_per_s,u0, v0

      double precision u0_comm,v0_comm, revs_comm
      common /comm_velocity/ u0_comm,v0_comm,revs_comm

      revs_per_s = revs_comm
      u0 = u0_comm
      v0 = v0_comm

      psi = -revs_per_s*(u0*xc - v0*yc)

      end
