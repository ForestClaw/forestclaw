      subroutine  rpn2_cons_update_manifold(meqn,maux, idir, q,aux,flux)

      implicit none

      integer meqn,maux,idir
      double precision q(meqn), aux(maux), flux(meqn)
      double precision u,v,urot

c        # Cell-centered velocities         
         u = aux(1)
         v = aux(2)

         urot = aux(6+2*idir)*u + aux(7+2*idir)*v
         flux(1) = urot*q(1)

      end
