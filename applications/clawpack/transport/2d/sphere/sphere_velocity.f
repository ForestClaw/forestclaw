c     # ------------------------------------------------------------
c     # Prescribes velocity fields for the unit sphere.
c     # 
c     # Assumes all components are given in coordinates relative to 
c     # the standard basis (1,0) and (0,1). 
c     # ------------------------------------------------------------
   
      subroutine velocity_components(x,y,t, u, vcart,flag)
      implicit none

      double precision x, y, t, u(2),vcart(3)
      double precision derivs(4)
      integer flag

      call velocity_derivs(x,y,t, u,vcart, derivs, flag)

      end

      subroutine velocity_derivs(x,y,t, u, vcart, derivs, flag)
      implicit none

      double precision x, y, t, u(2)
      double precision derivs(4)
      integer flag

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example        

      double precision omega(3)
      common /rotation_comm/ omega

      double precision kappa, period
      common /wind_comm/ kappa, period      

      double precision u1x, u1y, u2x, u2y
      double precision phi, theta

      double precision tp, lp
      double precision t1(3), t2(3), xp, yp, zp, w
      double precision rv(3), vcart(3)


c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      call map_comp2spherical(x,y,phi,theta)

      u1x = 0
      u1y = 0
      u2x = 0
      u2y = 0

      flag = 0

      if (example .eq. 0) then
          flag = 1
          call mapc2m_spherical(x,y,xp,yp,zp)
          rv(1) = xp
          rv(2) = yp
          rv(3) = zp
          call map_cross(omega,rv, vcart,w)

c         # Define (du(1)/dx, du(2)/dy, du(3)/dz)          
          derivs(1) = 0
          derivs(2) = 0
          derivs(3) = 0
      elseif (example .eq. 1) then
          flag = 0
          tp = t/period         
          lp = phi - pi2*tp

          u(1) = 10/period*sin(lp)**2*sin(2*theta)*cos(pi*tp) +
     &                    2*pi*cos(theta)/period
          u(2) = 10/period*sin(2*lp)*cos(theta)*cos(pi*tp)
          
c         u1x = -2*pi*s*cos(pi*x)*sin(pi*x)
c         u2y =  2*pi*s*sin(pi*y)*cos(pi*y)

      elseif (example .eq. 2) then
          flag = 0
          u(1) = -5/period*sin(lp/2.d0)**2*sin(2*theta)*
     &                    cos(theta)**2*cos(pi*tp) +
     &                    2*pi*cos(theta)/period
          u(2) = 5/(2*period)*sin(lp)*cos(theta)**3*cos(pi*tp)
      else
         write(6,'(A,A)') 'sphere_velocity : ',
     &              'No valid example provided'
         stop
      endif

      if (flag .eq. 0) then
          derivs(1) = u1x
          derivs(2) = u1y
          derivs(3) = u2x
          derivs(4) = u2y
      endif          

      end


      subroutine velocity_components_cart(x,y,t,vcart)
      implicit none

      double precision x,y,t, vcart(3)
      double precision u(2), t1(3), t2(3)
      integer flag, k


      call velocity_components(x,y,t, u,vcart,flag)

      if (flag .eq. 0) then
c         # Velocity components are given in Cartesian components
          call map_covariant_basis(x, y, t1,t2)

          do k = 1,3
              vcart(k) = u(1)*t1(k) + u(2)*t2(k)
          enddo
      endif

      end


c     # ------------------------------------------------------------
c     # Called from qexact
c     # ------------------------------------------------------------
      subroutine velocity_components_spherical(x,y,t,u)
      implicit none

      double precision x,y,t, u(2)
      double precision vcart(3), t1(3), t2(3)
      double precision t1n2, t2n2, map_dot
      integer flag

      call velocity_components(x,y,t, u,vcart,flag)

      if (flag .eq. 1) then
c         # Velocity components are given in Cartesian components
          call map_covariant_basis(x, y, t1,t2)
          t1n2 = map_dot(t1,t1)
          t2n2 = map_dot(t2,t2)
          u(1) = map_dot(vcart,t1)/t1n2
          u(2) = map_dot(vcart,t2)/t2n2
      endif

      end


c     # ------------------------------------------------------------
c     # Public interface (called from setaux)
c     # ------------------------------------------------------------
      subroutine sphere_center_velocity(x,y,t,vcart)
      implicit none

      double precision x,y,t, vcart(3)

      call velocity_components_cart(x,y,t,vcart)

      end




