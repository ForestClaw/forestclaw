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

      integer example
      common /example_comm/ example  

      double precision psi

      pi2 = 2*pi      

c     # Velocity is described in terms of a stream function.  This 
c     # approach depends on the mapping. 

      if (example .eq. 0) then

c         # Rigid body rotation
          if (mapping .eq. 0 .or. mapping .eq. 2) then
              psi = (pi2*revs_per_s)*
     &                alpha*(pi2*yc1 + alpha*sin(pi2*yc1))
          elseif (mapping .eq. 1) then
              psi = (pi2*revs_per_s)*alpha*
     &                  (pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
          endif
      else

          write(6,'(A,A)') '(psi.f) Non-rigid body rotation not ',
     &               'defined in terms of a streamfunction'
          stop
      endif
      torus_psi = psi

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

      integer example
      common /example_comm/ example  


      pi2 = 2*pi

      if (example .eq. 0) then
          if (mapping .eq. 0) then
c             # psi = (pi2*revs_per_s)*alpha*(pi2*yc1 + alpha*sin(pi2*yc1))
              psi_xi = 0
              psi_eta = (pi2)**2*revs_per_s*alpha*(1 + 
     &                    alpha*cos(pi2*yc1))
          elseif (mapping .eq. 1) then 
c             # psi = (pi2*revs_per_s)*alpha*(pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
              psi_xi = (pi2)**2*revs_per_s*alpha*(1 + 
     &           alpha*cos(pi2*(xc1+yc1)))
              psi_eta = (pi2)**2*revs_per_s*alpha*(1 + 
     &           alpha*cos(pi2*(xc1+yc1)))
    
          endif 
      else
         write(6,'(A,A)') '(psi.f) : Velocity not defined in terms', 
     &                    'of a stream function'
         stop
      endif


      end

      subroutine torus_velocity_components(xc1,yc1,u1,u2,u11,u22)
      implicit none

      double precision xc1, yc1, u1, u2, u11, u22
      double precision t1(3), t2(3), s1(3), s2(3)
      double precision m11, m12, m21, m22, m13, m23, v1, v2
      double precision minv11, minv12, minv21, minv22, det
      double precision minv13, minv23
      double precision nvec1(3), nvec2(3), w1 , w2
      double precision torus_dot

      double precision tau1(3), tau2(3), d1(3), d2(3)
      double precision xp,yp,zp,xc2, yc2
      integer blockno_dummy

      integer*8 cont, get_context

      double precision pi
      common /compi/ pi

      integer example
      common /example_comm/ example  

      integer mapping
      common /mapping_comm/ mapping

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision alpha
      common /torus_comm/ alpha

      double precision s
      integer k

      cont = get_context()

      if (example .eq. 0) then
c         # Rigid body rotation 
          if (mapping .eq. 0) then
c             # Counter clockwise?             
              u1 = revs_per_s
              u2 = 0
              u11 = 0
              u22 = 0
          elseif (mapping .eq. 1) then

c             # First, invert map to get (xc,yc) for this point in 
c             # the regular torus map.              
c             # Do not map (xc1,yc1) to brick, since mapping was done above
              blockno_dummy = -1  
              call fclaw2d_map_c2m(cont,blockno_dummy,xc1,yc1,xp,yp,zp)
              call mapc2m_torus_invert(xp,yp,zp,xc2,yc2,alpha)

              mapping = 0  !! Change value of common block variable
              call torus_contravariant_basis(xc2,yc2,t1,t2)
              mapping = 1  !! Change value back
              call torus_covariant_basis(xc1,yc1,s1,s2)   

              m11 = torus_dot(s1,t1)
              m12 = torus_dot(s1,t2)
              m21 = torus_dot(s2,t1)
              m22 = torus_dot(s2,t2)

              det = m11*m22 - m12*m21
              minv11 = m22/det
              minv12 = -m12/det
              minv21 = -m21/det
              minv22 = m11/det

              v1 = revs_per_s   !! from above
              v2 = 0
              u1 = minv11*v1 + minv21*v2
              u2 = minv12*v1 + minv22*v2
              u11 = 0
              u22 = 0
c             # Twisted torus
          endif

      elseif (example .eq. 1) then
c         # Some other velocity field (div u not necessarily zero)
          if (mapping .eq. 0) then
c             # Velocity field with divergence        
              s = sqrt(2.d0)
              u1 = s*cos(8*pi*xc1)
              u2 = s*sin(8*pi*yc1)   
    
              u11 = -8*pi*s*sin(8*pi*xc1)
              u22 = 8*pi*s*cos(8*pi*xc1)          
          elseif (mapping .eq. 1) then
c             #  velocity components defined in terms of twisted torus           
          endif   
      endif

      end




