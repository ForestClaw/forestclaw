c     # ------------------------------------------------------------
c     # Prescribes velocity fields for the unit sphere.
c     # 
c     # ------------------------------------------------------------
   
      subroutine cylinder_velocity_derivs(x,y,t, u,vcart,derivs,
     &                                     cart_flag)
      implicit none

      double precision x, y, u(2), vcart(3), t
      double precision derivs(4)
      integer cart_flag

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      double precision revs_per_s, cart_speed
      common /stream_comm/ revs_per_s, cart_speed

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /cylinder_comm/ alpha, beta, theta_range, phi_range

      double precision s, a, xc1, yc1
      double precision theta, thetax, thetay
      double precision phi, phix, phiy

      call map_comp2cylinder_derivs(x,y,theta,phi,thetax,thetay,
     &                             phix, phiy)

      if (example .eq. 0) then
          cart_flag = 0
          u(1) = revs_per_s
          u(2) = 0
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
          derivs(4) = 0
      elseif (example .eq. 1) then
          cart_flag = 0
          u(1) = 0
          u(2) = revs_per_s
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
          derivs(4) = 0
      elseif (example .eq. 2) then        
          cart_flag = 1
          vcart(1)  = cart_speed
          vcart(2) = 0
          vcart(3) = 0
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
      elseif (example .eq. 3) then
c         # This flow is divergent (not divergence-free)          
          cart_flag = 0
          u(1) = 1
          u(2) = 1
c         # uderivs = [u1x u1y; u2x u2y]          
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0 
          derivs(4) = 0
      endif


      end



c     # ----------------------------------------------------------------
c     #                       Public interface
c     # ----------------------------------------------------------------


c     # ---------------------------------------------
c     # Called from setaux
c     # 
c     #    -- used to compute velocity at faces
c     # ---------------------------------------------
      subroutine velocity_components_cart(x,y,t,vcart)
      implicit none

      double precision x,y,t, u(2), vcart(3), uderivs(4)
      double precision t1(3), t2(3)
      integer cart_flag, k


      call cylinder_velocity_derivs(x,y,t, u,vcart,uderivs,cart_flag)

      if (cart_flag .eq. 0) then
c         # Velocity components are given in spherical components
c         # and must be converted to Cartesian
          call map_covariant_basis(x, y, t1,t2)

          do k = 1,3
              vcart(k) = u(1)*t1(k) + u(2)*t2(k)
          enddo
      endif

      end


c     # ------------------------------------------------------------
c     # Called from map_divergence
c     # 
c     #    -- Needed to define ODE system to get exact solution
c     # ------------------------------------------------------------
      subroutine velocity_derivs(x,y,t, u, vcart, derivs, cart_flag)
      implicit none

      double precision x,y,t, u(2), vcart(3), derivs(4)
      double precision t1(3), t2(3), t1n2, t2n2, map_dot
      double precision t1inv(3), t2inv(3)
      integer cart_flag

      call cylinder_velocity_derivs(x,y,t, u,vcart,derivs,cart_flag)

      if (cart_flag .eq. 1) then
c         # Velocity components are given in Cartesian components
          call map_contravariant_basis(x, y, t1inv,t2inv)
          u(1) = map_dot(vcart,t1inv)
          u(2) = map_dot(vcart,t2inv)

c         # Need to convert derivatives to derivatives with respect
c         # to computational coordinates.          

      endif

      end

c     # ------------------------------------------------------------
c     # Called from qexact
c     # 
c     #  -- components relative to basis are needed.
c     # ------------------------------------------------------------
      subroutine velocity_components(x,y,t,u)
      implicit none

      double precision x,y,t, u(2), vcart(3), derivs(4)
      integer cart_flag

      call velocity_derivs(x,y,t, u,vcart,derivs,cart_flag)

      end



