c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=torus; 1=twisted torus)
c     # ------------------------------------------------------------
      double precision function torus_psi(xc1,yc1)
      implicit none

      double precision xc1, yc1

      double precision pi, pi2
      common /compi/ pi

      double precision alpha
      common /torus_comm/ alpha

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer mapping
      common /mapping_comm/ mapping

      double precision psi

      pi2 = 2*pi

c     # Velocity is described in terms of (xi, eta)=(xc1,yc1) coordinates
      if (mapping .eq. 0 .or. mapping .eq. 2) then
c         # Rigid body rotation
          psi = (pi2*revs_per_s)*alpha*(pi2*yc1 + alpha*sin(pi2*yc1))
      elseif (mapping .eq. 1) then
c         # Twisted torus stream function (to be used with usual torus map)
          psi = (pi2*revs_per_s)*alpha*
     &                (pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
      endif
      torus_psi = psi

      end

      subroutine torus_velocity_components(xc1,yc1,u1,u2,u11,u22)
      implicit none

      double precision xc1, yc1, u1, u2, u11, u22

      double precision pi
      common /compi/ pi

      integer example
      common /example_comm/ example  

      integer mapping
      common /mapping_comm/ mapping

      double precision s

      if (mapping .eq. 0) then
          if (example .eq. 0) then
c             # Rigid body rotation (assuming regular torus)
              u1 = 1.d0
              u2 = 0
              u11 = 0
              u22 = 0
          elseif (example .eq. 1) then
c             # Velocity field with divergence        
              s = sqrt(2.d0)
              u1 = s*cos(8*pi*xc1)
              u2 = s*sin(8*pi*yc1)   
    
              u11 = -8*pi*s*sin(8*pi*xc1)
              u22 = 8*pi*s*cos(8*pi*xc1)
          endif
      elseif (mapping .eq. 1) then
c         #  velocity components defined in terms of twisted torus           
      endif   

      end


      subroutine torus_psi_derivs(xc1,yc1,psi_xi, psi_eta)
      implicit none

      double precision xc1, yc1
      double precision psi_xi, psi_eta

      double precision pi, pi2
      common /compi/ pi

      double precision alpha
      common /torus_comm/ alpha

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer mapping
      common /mapping_comm/ mapping

      pi2 = 2*pi

      if (mapping .eq. 0) then
c         psi = (pi2*revs_per_s)*alpha*(pi2*yc1 + alpha*sin(pi2*yc1))
          psi_xi = 0
          psi_eta = (pi2)**2*revs_per_s*alpha*(1 + alpha*cos(pi2*yc1))
      elseif (mapping .eq. 1) then 
c         psi = (pi2*revs_per_s)*alpha*(pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))

          psi_xi = (pi2)**2*revs_per_s*alpha*(1 + 
     &           alpha*cos(pi2*(xc1+yc1)))
          psi_eta = (pi2)**2*revs_per_s*alpha*(1 + 
     &           alpha*cos(pi2*(xc1+yc1)))

      endif 


      end


