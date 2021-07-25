      subroutine  rpn2cons_update_manifold(meqn,maux, idir, iface, q,
     &                                     auxvec_center,
     &                                     auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, k, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision u,v,urot

c     # Cell-centered velocities         
      u = auxvec_center(1)
      v = auxvec_center(2)

c     # x-face normal : (6,7)
c     # y-face normal : (8,9)       
      if (idir .eq. 0) then
          urot = auxvec_edge(6)*u + auxvec_edge(7)*v
      else
          urot = auxvec_edge(8)*u + auxvec_edge(9)*v
      endif

c     #  f(q) = (n dot u)*q
      flux(1) = urot*q(1)

      end
