c     # ------------------------------------------------------------
c     # Compute cell-centered velocity fields
c     # 
c     # Define cell-centered velocity field using basis functions
c     # 
c     #      u = u1*t1 + u2*t2   (div u might not be zero)
c     # 
c     # NOTE: All arguments to routines here should be in computational 
c     # coordinates (x,y) in [0,1]x[0,1].  Mapping from brick domains
c     # to [0,1]x[0,1] should be done from calling routines, as should 
c     # any mappings that convert from an orthogonal coordinate system
c     # to a non-orthogonal system.
c     # ------------------------------------------------------------



      subroutine annulus_velocity_derivs(x,y,t, u,vcart,derivs,flag)
      implicit none

      double precision x, y, t, u(2), vcart(3), derivs(4)
      integer flag

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s, cart_speed, amplitude, freq
      common /velocity_comm/ revs_per_s, cart_speed, amplitude, freq

      integer example
      common /example_comm/ example  

      double precision xp,yp,zp, ravg, xc, d, tfinal, A
      double precision r, theta, w, nc

c     # Set non-zeros derivs only
      if (example .eq. 0) then
         flag = 0
c        # Rigid body rotation        
         u(1) = revs_per_s
         u(2) = 0
         derivs(1) = 0
         derivs(2) = 0
         derivs(3) = 0
         derivs(4) = 0
      else
          flag = 1
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
          if (example .eq. 1) then
              vcart(1) = cart_speed
              vcart(2) = 0
          elseif (example .eq. 2) then
              vcart(1) = 0
              vcart(2) = cart_speed
          elseif (example .eq. 3) then
             call map_comp2annulus(x,y,theta,r)
             w = 0.5
             vcart(1) = -(1-w)*pi2*r*sin(theta)
             vcart(2) = (1+w)*pi2*r*cos(theta)
             nc = sqrt(vcart(1)**2 + vcart(2)**2)

             vcart(1) = -vcart(1)/nc
             vcart(2) = -vcart(2)/nc
         elseif (example .eq. 4) then
              A = amplitude
              tfinal = 0.25
              vcart(1) = cart_speed
              vcart(2) = pi2*A*cos(freq*pi2*t/tfinal)/tfinal;
          elseif (example .eq. 5) then
              A = amplitude
              tfinal = 0.25
              vcart(1) = cart_speed*pi*sin(pi*t/tfinal)/2.d0
              vcart(2) = 0
          endif
          vcart(3) = 0
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


      call annulus_velocity_derivs(x,y,t, u,vcart,uderivs,flag)

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

      call annulus_velocity_derivs(x,y,t, u,vcart,derivs,flag)

      if (flag .eq. 1) then
c         # Velocity components are given in Cartesian components
          call map_covariant_basis(x, y, t1,t2)
          t1n2 = map_dot(t1,t1)
          t2n2 = map_dot(t2,t2)
          u(1) = map_dot(vcart,t1)/t1n2
          u(2) = map_dot(vcart,t2)/t2n2

c          call map_contravariant_basis(x, y, t1inv,t2inv)
c          # Convert Cartesian derivatives to u(1),u(2) derivatives
c          # 
c          #   grad v1 = (dv1/dx)*t1c + (dv1/dy)*t2c + (dv1/dz)*t3c
c          #   grad v2 = (dv2/dx)*t1c + (dv2/dy)*t2c + (dv2/dz)*t3c
c          # .....
c          # Don't think this is needed ...

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




