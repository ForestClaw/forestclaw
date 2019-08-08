c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=swirl; 1=twisted swirl)
c     # ------------------------------------------------------------


      subroutine velocity_components(blockno,x,y,u)
      implicit none

      double precision x, y, u(2)
      integer blockno
      double precision uderivs(4)

      double precision s

      call velocity_derivs(blockno, x,y,u,uderivs)

      end

      subroutine velocity_derivs(blockno,x,y,u,uderivs)
      implicit none

      double precision x, y, u(2)
      integer blockno
      double precision uderivs(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      double precision velocity(2)
      common /velocity_comm/ velocity      

      double precision s, pim, u1x, u1y, u2x, u2y
      integer k


c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      u1x = 0
      u1y = 0
      u2x = 0
      u2y = 0

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
         write(6,'(A,A)') 'clawpack46_setaux : ',
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
      subroutine square_center_velocity(blockno, x,y,vel)
      implicit none

      double precision x,y,vel(3)
      integer blockno

      double precision t1(3), t2(3)
      double precision t1inv(3), t2inv(3)
      double precision nvec(3), gradpsi(3), sv

      double precision p, px, py
      double precision u(2), uderivs(2)

      integer k

c     # Vector field defined as u1*tau1 + u2*tau2    

c      call map_covariant_basis(blockno, x, y, t1,t2)
      call velocity_components(blockno, x,y,u)


c     # Velocities are all given in terms of Cartesian components   
      do k = 1,2
        vel(k) = u(k)
      enddo
      vel(3) = 0
        
      end




