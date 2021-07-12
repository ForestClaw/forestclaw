c     # ------------------------------------------------------------
c     # Prescribes velocity fields for the unit square ([0.1]x[0,1]).
c     # 
c     # Assumes all components are given in coordinates relative to 
c     # the standard basis (1,0) and (0,1). 
c     # ------------------------------------------------------------


      subroutine square_velocity_derivs(x,y,t,u,vcart,derivs,flag)
      implicit none

      double precision x, y, t, u(2), vcart(3)
      double precision derivs(4)
      integer flag

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      double precision velocity(2)
      common /velocity_comm/ velocity      

      double precision s, u1x, u1y, u2x, u2y


c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      u1x = 0
      u1y = 0
      u2x = 0
      u2y = 0

c     # Basis map is the unit 
      flag = 1      

c     # Set non-zeros derivs only
      s = sqrt(2.d0)
      if (example .eq. 0) then
         u(1) = velocity(1)
         u(2) = velocity(2)
      elseif (example .eq. 1) then
c        # No sonic points, i.e. velocity field > 0
         u(1) = s*(cos(pi*x)**2 + 0.5d0)         
         u(2) = s*(sin(pi*y)**2 + 0.5d0)
         u1x = -2*pi*s*cos(pi*x)*sin(pi*x)
         u2y =  2*pi*s*sin(pi*y)*cos(pi*y)
      elseif (example .eq. 2) then
c        # Velocity field crosses 0
         u(1) = s*(cos(pi*x)**2 - 0.5d0)
         u(2) = s*(sin(pi*y)**2 - 0.5d0)
         u1x = -2*pi*s*cos(pi*x)*sin(pi*x)
         u2y =  2*pi*s*sin(pi*y)*cos(pi*y)
      else
         write(6,'(A,A)') 'square_velocity : ',
     &              'No valid example provided'
         stop
      endif

      vcart(1) = u(1)
      vcart(2) = u(2)
      vcart(3) = 0

      derivs(1) = u1x
      derivs(2) = u1y
      derivs(3) = u2x
      derivs(4) = u2y

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

      double precision x,y,t, u(2), vcart(3), derivs(4)
      integer flag

      call square_velocity_derivs(x,y,t, u,vcart,derivs,flag)

c      if (flag .eq. 0) then
cc         # Velocity components are given in spherical components
cc         # and must be converted to Cartesian
c          call map_covariant_basis(x, y, t1,t2)
c
c          do k = 1,3
c              vcart(k) = u(1)*t1(k) + u(2)*t2(k)
c          enddo
c      endif

      end


c     # ------------------------------------------------------------
c     # Called from map_divergence
c     # 
c     #    -- Needed to define ODE system to get exact solution
c     # ------------------------------------------------------------
      subroutine velocity_derivs(x,y,t, u, vcart, derivs, flag)
      implicit none

      double precision x,y,t, u(2), vcart(3), derivs(4)
      integer flag

      call square_velocity_derivs(x,y,t, u,vcart,derivs,flag)

c     # Not needed for Cartesian grid
c      if (flag .eq. 1) then
cc         # Velocity components are given in Cartesian components
cc         # Derivatives are automatically given in terms of basis.
c          call map_covariant_basis(x, y, t1,t2)
c          t1n2 = map_dot(t1,t1)
c          t2n2 = map_dot(t2,t2)
c          u(1) = map_dot(vcart,t1)/t1n2
c          u(2) = map_dot(vcart,t2)/t2n2
c      endif

      end

c     # ------------------------------------------------------------
c     # Called from qexact
c     # 
c     #  -- components relative to basis are needed.
c     # ------------------------------------------------------------
      subroutine user_velocity_components_cart(x,y,t,vcart)
      implicit none

      double precision x,y,t, u(2), vcart(3), derivs(4)
      integer flag

      call velocity_derivs(x,y,t, u,vcart,derivs,flag)

      end

      subroutine user_map2comp(blockno,xc,yc,xp,yp,zp,xc1,yc1)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp,xc1,yc1

      integer*8 cont, get_context

      double precision zc1

      cont = get_context()

      !! This doesn't work with the five patch mapping
      !! call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

      xc1 = xp
      yc1 = yp
      zc1 = 0


      end



