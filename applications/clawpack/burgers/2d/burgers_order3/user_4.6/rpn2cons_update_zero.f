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
