      subroutine  rpn2cons_update(meqn,maux, idir, iface, q,
     &                             auxvec_center,
     &                             auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, m, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision urot

      double precision ubar, vbar
      common /comrp/ ubar,vbar


      if (idir .eq. 0) then
         urot = ubar
      else
         urot = vbar
      endif

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

