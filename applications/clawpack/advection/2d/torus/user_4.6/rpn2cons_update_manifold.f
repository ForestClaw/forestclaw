      subroutine  rpn2_cons_update_manifold(meqn,maux, idir, iface, q,
     &                                     auxvec_center,
     &                                     auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, m, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)
      double precision u,v,w,urot, nv(3)
      integer k

cc     # Cell-centered velocities         
c      u = auxvec_edge(4)
c      v = auxvec_edge(5)
c      w = auxvec_edge(6)

c     # x-face normal : (9-11)
c     # y-face normal : (12-14)       
c      do k = 1,3
c          if (idir .eq. 0) then
c              nv(k) = auxvec_edge(9+k-1)
c          else
c              nv(k) = auxvec_edge(12+k-1)
c          endif
c      enddo

c     # Get flux in direction normal to idir-face
c      urot = u*nv(1) + v*nv(2) + w*nv(3)

      if (iface .eq. 0) then
          urot = auxvec_center(2+2*idir)  !! Left edge or bottom edge
      else
          urot = auxvec_center(3+2*idir)
      endif

c     #  f(q) = (n dot u)*q
      do m = 1,meqn
c         # Don't multiply by edgelength (scaling is done elsewhere)
          flux(m) = urot*q(m)
      enddo

      end

      subroutine  rpn2_cons_update_zero(meqn,maux, idir, q,
     &                                  auxvec_center,
     &                                  auxvec_edge,flux)

      implicit none

      integer meqn,maux,idir, m
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

c     #  f(q) = (n dot u)*q
      do m = 1,meqn
c         # Don't multiply by edgelength (scaling is done elsewhere)
          flux(m) = 0
      enddo

      end

