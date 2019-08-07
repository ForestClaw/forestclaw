      subroutine  rpn2cons_update(meqn,maux, idir, iface, q,
     &                             auxvec_center,
     &                             auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, m, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision fp

c     #  f(q) = (n dot u)*q

      fp = auxvec_center(2)
      do m = 1,meqn
c         # Don't multiply by edgelength (scaling is done elsewhere)
          flux(m) = fp
      enddo

      end

