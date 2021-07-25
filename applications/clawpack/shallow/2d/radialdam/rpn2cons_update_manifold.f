c     # This routine rotates cell-centered data into frame specified
c     # by vectors at the edge.  The flux is then computed from this 
c     # rotated data.   
c     # 
c     # All scaling and potential sign changes are handled by the 
c     # calling routine.
c
c     # Note that edge tangents must be normalized.

      subroutine  rpn2_cons_update_manifold(meqn,maux, idir, iface, q,
     &                                     auxvec_center,
     &                                     auxvec_edge,flux)
      implicit none

      integer meqn,maux,idir, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

      double precision grav
      common /cparam/  grav

      double precision enx, eny, enz
      double precision etx, ety, etz
      double precision hun, hut, h, un, ut, gamma, f(4)
      integer ioff

c     #  f1 = (hu; hu^2 + 0.5*g*h^2; huv)
c     #  f2 = (hv; huv; hv^2 + 0.5*gh^2)

      if (idir .eq. 0) then
          ioff = 1
      else
          ioff = 7
      endif

      enx =   auxvec_edge(ioff+1)
      eny =   auxvec_edge(ioff+2)
      enz =   auxvec_edge(ioff+3)      
      etx =   auxvec_edge(ioff+4)
      ety =   auxvec_edge(ioff+5)
      etz =   auxvec_edge(ioff+6)

c     !! Normalize the edge lengths
      gamma = dsqrt(etx**2 + ety**2 + etz**2)
      etx =   etx / gamma
      ety =   ety / gamma
      etz =   etz / gamma
      
      hun = enx*q(2)   + eny*q(3)   + enz*q(4)
      hut = etx*q(2)   + ety*q(3)   + etz*q(4)

      h = q(1)
      un = hun/h
      ut = hut/h

      f(1) = hun
      f(2) = hun**2/h + 0.5*grav*h**2
      f(3) = un*hut
      f(4)  = un*q(4)

      flux(1) = hun
      flux(2) = enx*f(2) + etx*f(3)
      flux(3) = eny*f(2) + ety*f(3)
      flux(4) = enz*f(2) + etz*f(3)

      end


