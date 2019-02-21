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

      subroutine torus_basic_map(xc1,yc1,u1,u2,coderiv1,
     &                           coderiv2,tau1,tau2)

      implicit none
      double precision xc1, yc1, u1, u2, tau1(3), tau2(3)
      double precision coderiv1(2), coderiv2(2)
      double precision u11, u22, u12, u21
      double precision torus_dot

      double precision pi
      common /compi/ pi

      integer example
      common /example_comm/ example  

      integer mapping, map_save
      common /mapping_comm/ mapping

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision alpha
      common /torus_comm/ alpha

      double precision s
      integer k

      map_save = mapping
      mapping = 0
      if (example .eq. 0) then
c         # Rigid body rotation 
c         # Counter clockwise?             
          u1 = revs_per_s
          u2 = 0

c             # Compute covariant derivatives
          u11 = 0
          u21 = 0
          call torus_covariant_derivative(xc1,yc1,1,u1,u2,
     &                                    u11,u21,coderiv1)

          u12 = 0
          u22 = 0
          call torus_covariant_derivative(xc1,yc1,2,u1,u2,
     &                                    u12,u22,coderiv2)

          call torus_covariant_basis(xc1,yc1,tau1,tau2)
          
      elseif (example .eq. 1) then
      endif

      mapping = map_save     
      end

      subroutine torus_velocity_components(xc1,yc1,u1,u2,coderiv1,
     &                                     coderiv2)
      implicit none

      double precision xc1, yc1, u1, u2
      double precision coderiv1(2), coderiv2(2)
      double precision coderiv01(2), coderiv02(2)
      double precision t1(3), t2(3), s1(3), s2(3)
      double precision m11, m12, m21, m22, v1, v2, det
      double precision minv11, minv22, minv12, minv21
      double precision u11, u22, u12, u21
      double precision torus_dot

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
              call  torus_basic_map(xc1,yc1,u1,u2,coderiv1,
     &                              coderiv2,t1,t2)
          elseif (mapping .eq. 1) then
c             # Twisted torus

c             # First, invert map to get (xc,yc) for this point in 
c             # the regular torus map.              
c             # Do not map (xc1,yc1) to brick, since mapping was done above
              blockno_dummy = -1  
              call fclaw2d_map_c2m(cont,blockno_dummy,xc1,yc1,xp,yp,zp)
              call mapc2m_torus_invert(xp,yp,zp,xc2,yc2,alpha)

              call  torus_basic_map(xc2,yc2,v1,v2,coderiv01,
     &                              coderiv02,t1,t2)

              call torus_contravariant_basis(xc1,yc1,s1,s2)   

c             # Tensor to handle change of basis
              m11 = torus_dot(t1,s1)
              m12 = torus_dot(t2,s1)
              m21 = torus_dot(t1,s2)
              m22 = torus_dot(t2,s2)

              u1 = m11*v1 + m12*v2
              u2 = m21*v1 + m22*v2

c             # Compute covariant derivatives.  Use mapping=0 basis
c             # and transform
              coderiv1(1) = m11*coderiv01(1) + m12*coderiv01(2)
              coderiv1(2) = m21*coderiv01(1) + m22*coderiv01(2)


              coderiv2(1) = m11*coderiv02(1) + m12*coderiv02(2)
              coderiv2(2) = m21*coderiv02(1) + m22*coderiv02(2)

c              write(6,100) coderiv01(1), coderiv01(2), 
c     &                 coderiv02(1), coderiv02(2)
c              write(6,100) coderiv1(1), coderiv1(2), 
c     &                 coderiv2(1), coderiv2(2)
c              write(6,*) ' '
100           format(4F24.16)              

          endif

      elseif (example .eq. 1) then
c         # Some other velocity field (div u not necessarily zero)
          if (mapping .eq. 0) then
c             # Velocity field with divergence        
              s = sqrt(2.d0)
              u1 = s*cos(8*pi*xc1)
              u2 = s*sin(8*pi*yc1)   


c             # Compute covariant derivatives
              u11 = -8*pi*s*sin(8*pi*xc1)
              u21 = 0
              call torus_covariant_derivative(xc1,yc1,1,u1,u2,
     &                                        u11,u21,coderiv1)

              u12 = 0
              u22 = 8*pi*s*cos(8*pi*xc1) 
              call torus_covariant_derivative(xc1,yc1,2,u1,u2,
     &                                        u12,u22,coderiv2)

              call torus_covariant_basis(xc1,yc1,t1,t2)
    
                       
          elseif (mapping .eq. 1) then
c             #  velocity components defined in terms of twisted torus           
          endif   
      endif

      end





