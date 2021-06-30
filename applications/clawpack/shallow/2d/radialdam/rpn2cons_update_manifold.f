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
      double precision hun, hut, h, un, ut, gamma
      integer ioff, mu, mv

c     #  f1 = (hu; hu^2 + 0.5*g*h^2; huv)
c     #  f2 = (hv; huv; hv^2 + 0.5*gh^2)

      if (idir .eq. 0) then
          ioff = 1
          mu = 2
          mv = 3
      else
          ioff = 7
          mu = 3
          mv = 2
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
      un = hun/q(1)
      ut = hut/q(1)

      flux(1) = hun
      flux(mu) = hun**2/h + 0.5*grav*h**2
      flux(mv) = hun*ut
      flux(4) = q(4)*un

      end

      subroutine  rpn2_cons_update(meqn,maux, idir, iface, q,
     &                             auxvec_center,
     &                             auxvec_edge,flux)
      implicit none

      integer meqn,maux,idir, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

      double precision grav
      common /cparam/  grav

      double precision hun, hut, h, un, ut

      integer mu, mv

c     #  f1 = (hu; hu^2 + 0.5*g*h^2; huv)
c     #  f2 = (hv; huv; hv^2 + 0.5*gh^2)

      if (idir .eq. 0) then
          mu = 2
          mv = 3
      else
          mu = 3
          mv = 2
      endif


      
      hun = q(mu)
      hut = q(mv)

      h = q(1)
      un = hun/h
      ut = hut/h

      flux(1) = hun
      flux(mu) = hun**2/h + 0.5*grav*h**2
      flux(mv) = hun*ut
      flux(4) = q(4)*un

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

