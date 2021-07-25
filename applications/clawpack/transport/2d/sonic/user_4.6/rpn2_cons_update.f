      subroutine  rpn2_cons_update(meqn,maux, idir, iface, 
     &                         q,auxvec_center,auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision u
      integer m

c     # Cell-centered velocities         
      u = auxvec_center(2 + iface)

      do m = 1,meqn
          flux(m) = u*q(m)
      end do

      end
