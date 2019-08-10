c     # ------------------------------------------------------------
c     # Prescribes velocity fields for the unit sphere.
c     # 
c     # Assumes all components are given in coordinates relative to 
c     # the standard basis (1,0) and (0,1). 
c     # ------------------------------------------------------------


      subroutine velocity_components(x,y,t, u)
      implicit none

      double precision x, y, t, u(2)
      double precision uderivs(4)

      call velocity_derivs(x,y,t, u,uderivs)

      end

      subroutine velocity_derivs(x,y,t, u,uderivs)
      implicit none

      double precision x, y, t, u(2)
      double precision uderivs(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example        

      double precision kappa, period
      common /wind_comm/ kappa, period      

      double precision kdiv2, kdiv4, tp, lp
      double precision u1x, u1y, u2x, u2y
      double precision phi, phix, phiy
      double precision theta, thetax, thetay
      integer k


c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      call map_comp2angles(x,y,phi,theta)

      u1x = 0
      u1y = 0
      u2x = 0
      u2y = 0

      kdiv2 = kappa/2.d0
      kdiv4 = kappa/4.d0
      tp = t/period

      lp = phi - pi2*tp

      if (example .eq. 0) then
c         u(1) = -kdiv2*sin(lp/2.d0)**2*sin(2*theta)*cos(pi*tp) +
c     &                2*pi*cos(theta)/period
c         u(2) = kdiv4*sin(lp)*(cos(theta)**3)*cos(pi*tp)

          u(1) = 1
          u(2) = 0

c         u1x = -2*pi*s*cos(pi*x)*sin(pi*x)
c         u2y =  2*pi*s*sin(pi*y)*cos(pi*y)
      else
         write(6,'(A,A)') 'sphere_velocity : ',
     &              'No valid example provided'
         stop
      endif

      uderivs(1) = u1x
      uderivs(2) = u1y
      uderivs(3) = u2x
      uderivs(4) = u2y

      end


c     # ------------------------------------------------------------
c     # Center : u = u1*tau1 + u2*tau2   (div u might not be zero)
c     # ------------------------------------------------------------
      subroutine sphere_center_velocity(x,y,t, vel)
      implicit none

      double precision x,y,t, vel(3)

      double precision t1(3), t2(3), u(2)

      integer k

c     # Vector field defined as u1*tau1 + u2*tau2    

      call map_covariant_basis(x, y, t1,t2)
      call velocity_components(x,y,t, u)

c     # Velocities are all given in terms of Cartesian components   
      do k = 1,3
        vel(k) = u(1)*t1(k) + u(2)*t2(k)
      enddo
        
      end




