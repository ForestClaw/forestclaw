      subroutine  rpn2_cons_update_manifold(meqn,maux, idir, iface, q,
     &                                     auxvec_center,
     &                                     auxvec_edge,flux)
      implicit none

      integer meqn,maux,idir, m, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision urot

c     # Get cell-centered velocity projected to face 
c     # 'iface' (in 0,1,2,3) 
      urot = auxvec_center(2+iface)

c     #  f(q) = (n dot u)*q
      do m = 1,meqn
c         # Don't multiply by edgelength (scaling is done elsewhere)
          flux(m) = urot*q(m)
      enddo

      end

      subroutine  rpn2_cons_update_zero(meqn,maux, idir, iface, q,
     &                                  auxvec_center,
     &                                  auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      integer m

c     #  f(q) = (n dot u)*q
      do m = 1,meqn
c         # No flux function available for equations in non-conservative form
          flux(m) = 0
      enddo

      end

