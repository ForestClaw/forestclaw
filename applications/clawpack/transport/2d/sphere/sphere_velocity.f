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
      double precision rv(3), vcart(3), th, t1n2, t2n2
      double precision map_dot, thetax, thetay, phix, phiy

      double precision lpx, lpy, thx, thy
      double precision cu1, fu1, fu1x, fu1y, gu1, gu1x, gu1y
      double precision cu2, fu2, fu2x, fu2y, gu2, gu2x, gu2y
      double precision hu1, hu1x, hu1y, hu2, hu2x, hu2y
      double precision gh, ghx, ghy
      double precision zf, zfx, zfy
      double precision uderivs_comp(4), ucomp(2)
      logical zonal_flow

      integer k

c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      call map_comp2spherical_derivs(x,y,theta,phi,
     &             thetax, thetay, phix, phiy)

      derivs(1) = 0
      derivs(2) = 0
      derivs(3) = 0
      derivs(4) = 0

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
      elseif (example .eq. 1 .or. example .eq. 2) then
c         # Divergence-free        
          flag = 0

c         # Lauritzen : (lambda,th)  <----> (theta,phi) (here)
c         # kappa = 10/period
          tp = t/period         
          lp = theta - pi2*tp
          th = phi  

          lpx = thetax
          lpy = thetay

          thx = phix
          thy = phiy

          if (example .eq. 1) then

c             # --------------------------------------------------
c             # From Nair and Lauritzen 2010 (Case 4; kappa = 2) 
c             # u = kappa*sin(lp)**2*sin(2*th)*cos(pi*t/Tfinal) +
c             #           2*pi*cos(th)/Tfinal
c             # v = kappa*sin(2*lp)*cos(th)*cos(pi*t/Tfinal)
c             # --------------------------------------------------

c             # u = U(1) component 
              cu1 = 10.0/period*cos(pi*tp)

              fu1 = sin(lp)**2
              fu1x = 2*sin(lp)*cos(lp)*lpx
              fu1y = 2*sin(lp)*cos(lp)*lpy

              gu1 = sin(2*th)
              gu1x = 2*cos(2*th)*thx
              gu1y = 2*cos(2*th)*thy

              ucomp(1) = cu1*fu1*gu1 + pi2*cos(th)/period
 
              u1x = cu1*(fu1*gu1x + fu1x*gu1) + 
     &                  -(pi2/period)*sin(th)*thx

              u1y = cu1*(fu1*gu1y + fu1y*gu1) + 
     &                  -(pi2/period)*sin(th)*thy


c             # v = U(2) component 
              cu2 = 10.0/period*cos(pi*tp)

              fu2 = sin(2*lp)
              fu2x = 2*cos(2*lp)*lpx
              fu2y = 2*cos(2*lp)*lpy

              gu2 = cos(th)
              gu2x = -thx*sin(th)
              gu2y = -thy*sin(th)

              ucomp(2) = cu2*fu2*gu2

              u2x = cu2*(fu2*gu2x + fu2x*gu2)
              u2y = cu2*(fu2*gu2y + fu2y*gu2)

          elseif (example .eq. 2) then
cc            # Velocity of spherical coordinates (from my original code)
c             u = -(kappa/2.d0)*sin(lp/2.d0)**2*sin(2*th)*cos(pi*t/Tfinal) +
c    &                2*pi*cth/Tfinal
c             v = (kappa/4.d0)*sin(lp)*(cth**3)*cos(pi*t/Tfinal)


c             # -----------------------------------------------------
c             # Nair and Lauritzen (2010) Case 3 (kappa = 2)   
c             # with background zonal flow         
c             #
c             u = -(kappa/2.d0)*sin(lp/2.d0)**2*sin(2*th)*cos(th)**2
c                    *cos(pi*t/Tfinal) + 2*pi*cos(th)/Tfinal
c             v = (kappa/4.d0)*sin(lp)*(cos(th)**3)*cos(pi*t/Tfinal)
c             # -----------------------------------------------------

c             !! From 2010 paper (no lp; no zonal flow)
              zonal_flow = .false.
              lp = theta              
  
c             # U(1) component 
              cu1 = -5.d0/period*cos(pi*tp)  

              fu1 = sin(lp/2.d0)**2
              fu1x = 2.d0*sin(lp/2)*cos(lp/2)*lpx/2.d0
              fu1y = 2.d0*sin(lp/2)*cos(lp/2)*lpy/2.d0

              gu1 = sin(2*th)
              gu1x = 2.d0*cos(2*th)*thx
              gu1y = 2.d0*cos(2*th)*thy

              hu1  = cos(th)**2
              hu1x = -2.d0*cos(th)*sin(th)*thx
              hu1y = -2.d0*cos(th)*sin(th)*thy


              if (zonal_flow) then
                  zf = (pi2/period)*cos(th)
                  zfx = -(pi2/period)*sin(th)*thx
                  zfy = 0
              else
                  zf = 0
                  zfx = 0
                  zfy = 0
              endif

              gh = gu1*hu1
              ghx = gu1*hu1x + gu1x*hu1
              ghy = gu1*hu1y + gu1y*hu1

              ucomp(1) = cu1*fu1*gu1*hu1 + zf

              u1x = cu1*(fu1*ghx + fu1x*gh) + zfx
              u1y = cu1*(fu1*ghy + fu1y*gh) + zfy


c             # U(2) component 
              cu2 = 5.d0/(2.d0*period)*cos(pi*tp)

              fu2 = sin(lp)
              fu2x = cos(lp)*lpx
              fu2y = cos(lp)*lpy

              gu2 = cos(th)**3
              gu2x = -3.d0*cos(th)**2*sin(th)*thx
              gu2y = -3.d0*cos(th)**2*sin(th)*thy

              ucomp(2) = cu2*fu2*gu2

              u2x = cu2*(fu2*gu2x + fu2x*gu2)
              u2y = cu2*(fu2*gu2y + fu2y*gu2)

          endif
c         # ------------------------------
c         # Store derivatives
c         # ------------------------------

          uderivs_comp(1) = u1x
          uderivs_comp(2) = u1y
          uderivs_comp(3) = u2x
          uderivs_comp(4) = u2y

c         # compute d(u1/norm(t1))/dx
c         # compute d(u1/norm(t1))/dy
c         # compute d(u2/norm(t2))/dx
c         # compute d(u2/norm(t2))/dy
          call  map_diff_normalized(x,y,ucomp, uderivs_comp, 
     &                                  derivs)              

c         # Normalize covariant vectors          
          call map_covariant_basis(x, y, t1,t2)
          t1n2 = map_dot(t1,t1)
          t2n2 = map_dot(t2,t2)
          u(1) = ucomp(1)/sqrt(t1n2)
          u(2) = ucomp(2)/sqrt(t2n2)
      else
         write(6,'(A,A)') 'sphere_velocity : ',
     &              'No valid example provided'
         stop
      endif

      end


      subroutine velocity_components_cart(x,y,t,vcart)
      implicit none

      double precision x,y,t, vcart(3)
      double precision u(2), t1(3), t2(3)
      integer flag, k


      call velocity_components(x,y,t, u,vcart,flag)

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




