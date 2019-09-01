c     # ------------------------------------------------------------
c     # Prescribes velocity fields for the unit sphere.
c     # 
c     # ------------------------------------------------------------
   
      subroutine torus_velocity_derivs(x,y,t, u,vcart,derivs,flag)
      implicit none

      double precision x, y, u(2), vcart(3), t
      double precision derivs(4)
      integer flag

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      double precision revs_per_s, cart_speed
      common /stream_comm/ revs_per_s, cart_speed

      double precision s, a, xc1, yc1
      double precision theta, thetax, thetay
      double precision phi, phix, phiy

      call map_comp2torus_derivs(x,y,theta,phi,thetax,thetay,
     &                             phix, phiy)

      if (example .eq. 0) then
          flag = 0
          u(1) = revs_per_s
          u(2) = 0
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
          derivs(4) = 0
      elseif (example .eq. 1) then
          flag = 0
          u(1) = 0
          u(2) = revs_per_s
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
          derivs(4) = 0
      elseif (example .eq. 2) then        
          flag = 1
          vcart(1)  = cart_speed
          vcart(2) = 0
          vcart(3) = 0
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
      elseif (example .eq. 3) then
          flag = 0
          s = sqrt(2.d0)
          a = 4
          u(1) = s*cos(a*theta)
          u(2) = s*sin(a*phi)   
c         # uderivs = [u1x u1y; u2x u2y]          
          derivs(1) = -s*a*sin(a*theta)*thetax
          derivs(2) = 0;
          derivs(3) = 0; 
          derivs(4) = s*a*cos(a*phi)*phiy
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
      integer flag, k


      call torus_velocity_derivs(x,y,t, u,vcart,uderivs,flag)

      if (flag .eq. 0) then
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
      subroutine velocity_derivs(x,y,t, u, vcart, derivs, flag)
      implicit none

      double precision x,y,t, u(2), vcart(3), derivs(4)
      double precision t1(3), t2(3), t1n2, t2n2, map_dot
      double precision t1inv(3), t2inv(3)
      integer flag

      call torus_velocity_derivs(x,y,t, u,vcart,derivs,flag)

      if (flag .eq. 1) then
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
      integer flag

      call velocity_derivs(x,y,t, u,vcart,derivs,flag)

      end



