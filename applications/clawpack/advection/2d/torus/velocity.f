c     # ------------------------------------------------------------
c     # Compute edge centered and cell-centered velocity fields
c     # 
c     # Edge : Compute velocity from a stream function.  
c     # This defines the average velocity at each edge
c     # 
c     #      u = curl \Psi  (div u = 0)
c     # 
c     # Center : Defines cell-centered velocities from a streamfunction
c     # 
c     #      u = n cross grad \Psi  (div u = 0)
c     # 
c     # or using basis functions
c     # 
c     #      u = u1*tau1 + u2*tau2   (div u might not be zero)
c     # 
c     # NOTE: All arguments to routines here should be in computational 
c     # coordinates (xi,eta) in [0,1]x[0,1].  Mapping from brick domains
c     # to [0,1]x[0,1] should be done from calling routines.
c     # ------------------------------------------------------------


c     # ------------------------------------------------------------
c     # Edge : u = \nabla cross \Psi - curl \Psi
c     # ------------------------------------------------------------
      subroutine torus_edge_velocity(xc1,yc1,xc2,yc2,ds,vn)
      implicit none

      double precision xc1,yc1,xc2,yc2, ds, vn, torus_psi

      vn = (torus_psi(xc1,yc1) - torus_psi(xc2,yc2))/ds

      end

c     # ------------------------------------------------------------
c     # Center : u = n cross grad \Psi   (div u = 0)
c     # 
c     # Center : u = u1*tau1 + u2*tau2   (div u might not be zero)
c     # ------------------------------------------------------------
      subroutine torus_center_velocity(xc1,yc1,vel)
      implicit none

      double precision xc1,yc1,vel(3)

      double precision tau1(3), tau2(3)
      double precision tau1inv(3), tau2inv(3)
      double precision nvec(3), gradpsi(3), sv
      logical use_stream

      double precision psi_xi, psi_eta
      double precision u1, u2, coderiv1(2), coderiv2(2)     

      double precision pi
      common /compi/ pi

      integer example
      common /example_comm/ example  

      integer k
    
      use_stream = .false.

      if (example .eq. 0 .and. use_stream) then
c         # Divergence free velocity field : u = n cross \Psi  
          call torus_psi_derivs(xc1,yc1,psi_xi,psi_eta)

          call torus_contravariant_basis(xc1,yc1, tau1inv,tau2inv)


          do k = 1,3
              gradpsi(k) = psi_xi*tau1inv(k) + psi_eta*tau2inv(k)
          end do

c         # Get surface normal
          call torus_covariant_basis(xc1,yc1,tau1,tau2)

c         Outward directed normal          
          call torus_cross(tau2,tau1,nvec,sv)

c         # Normalize surface normal
          do k = 1,3
              nvec(k) = nvec(k)/sv
          end do

c         # v = nvec x grad \Psi
          call torus_cross(nvec,gradpsi,vel,sv)
      else
c         # Vector field defined as u1*tau1 + u2*tau2    

          call torus_covariant_basis(xc1, yc1, tau1,tau2)
          call torus_velocity_components(xc1,yc1,u1,u2, coderiv1,
     &                                   coderiv2)

          do k = 1,3
              vel(k) = u1*tau1(k) + u2*tau2(k)
          enddo
      endif
        
      end


c     # ----------------------------------------------------------------
c     # Mapping functions 
c     # 
c     # Covariant basis vectors       : tau1 = T_xi, tau2 = T_eta
c     #
c     # Contravariant basis vectors   : tau1inv = a^{11}T_xi + a^{12}T_eta 
c     #                                 tau2inv = a^{21}T_xi + a^{22}T_eta 
c     # 
c     # Christoffel symbols           : g111, g112, g212, g222 (needed for 
c     #                                 computation of the divergence) 
c     #
c     # Mapping function depends on these global variables
c     #
c     #         mapping (0=torus; 1=twisted torus)
c     # ----------------------------------------------------------------

      subroutine torus_covariant_basis(xc1,yc1,tau1,tau2)
      implicit none

      double precision xc1, yc1
      double precision tau1(3), tau2(3)

      double precision R, Reta, Rxi


      double precision pi, pi2
      common /compi/ pi

      double precision alpha
      common /torus_comm/ alpha

      integer mapping
      common /mapping_comm/ mapping

      pi2 = 2*pi

      if (mapping .eq. 0) then
          R = 1 + alpha*cos(pi2*yc1)
          Reta = -pi2*alpha*sin(pi2*yc1)

c         # T_xi
          tau1(1) = -pi2*R*sin(pi2*xc1)
          tau1(2) = pi2*R*cos(pi2*xc1)
          tau1(3) = 0

c         # T_eta
          tau2(1) = Reta*cos(pi2*xc1)
          tau2(2) = Reta*sin(pi2*xc1)
          tau2(3) = pi2*alpha*cos(pi2*yc1)               

      elseif (mapping .eq. 1) then
c         # Coordinate normals          
          R    = 1 +  alpha*cos(pi2*(xc1 + yc1))
          Rxi  = -pi2*alpha*sin(pi2*(xc1 + yc1))
          Reta = -pi2*alpha*sin(pi2*(xc1 + yc1))

c         # T_xi
          tau1(1) = Rxi*cos(pi2*xc1) - pi2*R*sin(pi2*xc1)
          tau1(2) = Rxi*sin(pi2*xc1) + pi2*R*cos(pi2*xc1)
          tau1(3) = pi2*alpha*cos(pi2*(xc1 + yc1))

c         # T_eta
          tau2(1) = Reta*cos(pi2*xc1)
          tau2(2) = Reta*sin(pi2*xc1)
          tau2(3) = pi2*alpha*cos(pi2*(xc1+yc1))
      endif


      end


      subroutine torus_contravariant_basis(xc1,yc1, tau1inv, tau2inv)
      implicit none

      double precision xc1,yc1

      double precision tau1(3), tau2(3)
      double precision a11, a22, a12, a21, det
      double precision a11inv, a22inv, a12inv, a21inv
      double precision tau1inv(3), tau2inv(3)
      double precision torus_dot

      integer k

      call torus_covariant_basis(xc1,yc1, tau1,tau2)

c     # Compute grad psi(xi,eta) 
      a11 = torus_dot(tau1,tau1)
      a22 = torus_dot(tau2,tau2)
      a12 = torus_dot(tau1,tau2)
      a21 = a12

c     # Determinant
      det = a11*a22 - a12*a21

c     # Contravariant vectors
      a11inv = a22/det
      a22inv = a11/det
      a12inv = -a12/det
      a21inv = -a21/det     
      do k = 1,3
          tau1inv(k) = a11inv*tau1(k) + a12inv*tau2(k)
          tau2inv(k) = a21inv*tau1(k) + a22inv*tau2(k)
      end do

      end

cc     # Derivative in direction x^i = (xi, eta) 
      subroutine torus_covariant_derivative(xc1,yc1,j,u1,u2,
     &                                  u1j,u2j,coderiv)
      implicit none
 
      double precision xc1, yc1, u1,u2, u1j, u2j, coderiv(2)
      integer j

      double precision g111,g112,g211,g212
      double precision g121,g122,g221,g222
      double precision torus_christoffel_sym

      double precision u11, u12, u21, u22

c      call torus_velocity_components(xc1,yc1,u1,u2,u11,u22,u12,u21)

c     # Christoffel symbol gijk
c     # i - component  (u1 or u2)
c     # j - derivative (xi or eta)
c     # k - dummy variable
      if (j .eq. 1) then
          u11 = u1j
          u21 = u2j
          g111 = torus_christoffel_sym(xc1,yc1,1,1,1)
          g112 = torus_christoffel_sym(xc1,yc1,1,1,2)
          coderiv(1) = u11 + (u1*g111 + u2*g112)

          g211 = torus_christoffel_sym(xc1,yc1,2,1,1)
          g212 = torus_christoffel_sym(xc1,yc1,2,1,2)
          coderiv(2) = u21 + (u1*g211 + u2*g212)

c          write(6,100) coderiv(1), coderiv(2)
      else
          u12 = u1j
          u22 = u2j
          g121 = torus_christoffel_sym(xc1,yc1,1,2,1)
          g122 = torus_christoffel_sym(xc1,yc1,1,2,2)
          coderiv(1) = u12 + (u1*g121 + u2*g122)

          g221 = torus_christoffel_sym(xc1,yc1,2,2,1)
          g222 = torus_christoffel_sym(xc1,yc1,2,2,2)
          coderiv(2) = u22 + (u1*g221 + u2*g222)

c          write(6,100) coderiv(1), coderiv(2)
c          write(6,100) g121, g221
      endif

100   format(3F24.16)
      end


      double precision function torus_christoffel_sym(xc1,yc1,i,j,k) 
      implicit none

      double precision xc1, yc1
      integer i,j,k
      double precision gijk

      double precision T1(3), T2(3)
      double precision T11(3), T22(3), T12(3), T21(3)
      double precision gi(3), gjk(3)

      double precision R, R2, R22
      double precision torus_dot 

      double precision pi, pi2
      common /compi/ pi

      double precision alpha
      common /torus_comm/ alpha

      integer mapping
      common /mapping_comm/ mapping

      integer kk

      pi2 = 2*pi

      if (mapping .eq. 0) then
          R   = 1 + alpha*cos(pi2*yc1)
          R2  = -pi2*alpha*sin(pi2*yc1)
          R22 = -pi2**2*alpha*cos(pi2*yc1)


cc        # T_xi
c         T1(1) = -pi2*R*sin(pi2*xc1)
c         T1(2) = pi2*R*cos(pi2*xc1)
c         T1(3) = 0

cc        # T_eta
c         T2(1) = Reta*cos(pi2*xc1)
c         T2(2) = Reta*sin(pi2*xc1)
c         T2(3) = pi2*alpha*cos(pi2*yc1)               

          T11(1) = -pi2**2*R*cos(pi2*xc1)
          T11(2) = -pi2**2*R*sin(pi2*xc1)
          T11(3) = 0

          T12(1) = -pi2*R2*sin(pi2*xc1)
          T12(2) =  pi2*R2*cos(pi2*xc1)
          T12(3) = 0

          T21(1) = -pi2*R2*sin(pi2*xc1)
          T21(2) =  pi2*R2*cos(pi2*xc1)
          T21(3) = 0

          T22(1) = R22*cos(pi2*xc1)
          T22(2) = R22*sin(pi2*xc1)
          T22(3) = 0
      elseif (mapping .eq. 1) then
c         # To do ....        
      endif

      call torus_covariant_basis(xc1,yc1,T1,T2)

      do kk = 1,3
          if (i .eq. 1) then
              gi(kk) = T1(kk)
          else
              gi(kk) = T2(kk)
          endif

          if (j .ne. k) then
              if (j .eq. 1) then
                  gjk(kk) = T12(kk)
              else
                  gjk(kk) = T21(kk)
              endif
          else
              if (j .eq. 1) then
c                 # j == k == 1                  
                  gjk(kk) = T11(kk)
              else
c                 # j == k == 2                  
                  gjk(kk) = T22(kk)
              endif
          endif
      end do

      gijk = torus_dot(gi,gjk)

      torus_christoffel_sym = gijk

      end

      double precision function torus_divergence(xc1,yc1)
      implicit none

      double precision xc1,yc1

      double precision divu, u1, u2, coderiv1(2), coderiv2(2)

c      call torus_covariant_derivative(xc1,yc1,1,coderiv1)
c      call torus_covariant_derivative(xc1,yc1,2,coderiv2)

      call torus_velocity_components(xc1,yc1,u1,u2,coderiv1,
     &                                     coderiv2)
      

      torus_divergence = coderiv1(1) + coderiv2(2)

      end


c     # ----------------------------------------------------------------
c     # RHS functions for ODE solver DOPRI5  (original ode45)
c     #        -- divfree field
c     #        -- field with divergence
c     # ----------------------------------------------------------------

      subroutine torus_rhs_divfree(n,t,sigma,f,rpar,ipar)
      implicit none

      integer n, ipar
      double precision t, sigma(n), f(n), rpar


      integer m, i
      double precision xc1,yc1, q, u1,u2
      double precision psi_xi, psi_eta
      double precision tau1(3), tau2(3), t1xt2(3), w

      xc1 = sigma(1)
      yc1 = sigma(2)

      call torus_covariant_basis(xc1,yc1,tau1,tau2)
      call torus_psi_derivs(xc1,yc1,psi_xi,psi_eta)

c     # Compute u dot grad q in computational coordinates
      call torus_cross(tau1,tau2,t1xt2,w);

c     # Evolve contravariant components of velocity field          
      u1 = -psi_eta/w
      u2 = psi_xi/w

c     # This is a clockwise velocity field ? 
      f(1) = u1
      f(2) = u2

      end

      subroutine torus_rhs_nondivfree(n,t,sigma,f,rpar,ipar)
      implicit none

      integer n, ipar
      double precision t, sigma(n), f(n), rpar


      double precision xc1,yc1, q
      double precision u1,u2,u11,u22
      double precision divu, torus_divergence
      double precision coderiv1(2), coderiv2(2)

c     # Track evolution of these three quantities
      xc1 = sigma(1)
      yc1 = sigma(2)
      q = sigma(3)

      call torus_velocity_components(xc1,yc1,u1,u2,coderiv1,coderiv2)

      divu = coderiv1(1) + coderiv2(2)
      if (abs(divu) .gt. 1e-10) then
c          write(6,100) coderiv1(1), coderiv2(2), divu
      endif
100   format(3F24.16)

      f(1) = u1
      f(2) = u2
      f(3) = -divu*q   !! Non-conservative case

      end

c     # ----------------------------------------------------------------
c     # Utility functions
c     # ----------------------------------------------------------------

      subroutine torus_cross(u,v,uxv,w)
      implicit none

      double precision u(3),v(3),uxv(3),w, torus_dot
      integer k

      uxv(1) =   u(2)*v(3) - u(3)*v(2)
      uxv(2) = -(u(1)*v(3) - u(3)*v(1))
      uxv(3) =   u(1)*v(2) - u(2)*v(1)

      w =  torus_dot(uxv,uxv)
      w = sqrt(w)      

      end

      double precision function torus_dot(u,v)
      implicit none

      double precision u(3),v(3)

      torus_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      end





