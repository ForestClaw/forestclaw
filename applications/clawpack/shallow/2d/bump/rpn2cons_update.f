c     # This routine computes the flux of cell-centered 
c     # 
c     # All scaling and potential sign changes are handled by the 
c     # calling routine.
c

      subroutine  rpn2_cons_update(meqn,maux, idir, iface, q,
     &                             auxvec_center,
     &                             auxvec_edge,flux)
      implicit none

      integer meqn,maux,idir, iface
      double precision q(meqn), flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

      double precision grav
      common /cparam/  grav

      double precision h, hun, hut, un, ut

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
      flux(mv) = un*hut

      end


