      subroutine  rpn2_cons_update(meqn,maux, idir, q,aux,flux)

      implicit none

      integer meqn,maux,idir
      double precision q(meqn), aux(maux), flux(meqn)
      double precision u

      integer iface

c     # Cell-centered velocities         
      iface = idir+1
      u = aux(iface)

      flux(1) = u*q(1)

      end
