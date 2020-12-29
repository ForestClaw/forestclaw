      subroutine  rpn2cons_update_order2(meqn,maux, idir, iface, q,
     &                             auxvec_center, auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, m, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

c     #  f(q) = (n dot u)*q

      do m = 1,meqn
c         # Don't multiply by edgelength (scaling is done elsewhere)
          flux(m) = 0.5d0*q(1)**2
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

