c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=annulus; 1=twisted annulus)
c     # ------------------------------------------------------------



      subroutine annulus_velocity_components(x,y,t,u)
      implicit none

      double precision x, y, u(2), t

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s, cart_speed, amplitude, freq
      common /stream_comm/ revs_per_s, cart_speed, amplitude, freq

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer example
      common /example_comm/ example  

      double precision t1(3), t2(3)
      double precision t1_norm2, t2_norm2
      double precision vcart(3), annulus_dot
      double precision t1_dot_vcart, t2_dot_vcart
      double precision xp,yp,zp, ravg, xc, d, tfinal, A
      double precision r, w, th, nc


      call annulus_covariant_basis(x, y, t1,t2) 

      t1_norm2 = annulus_dot(t1,t1)
      t2_norm2 = annulus_dot(t2,t2)

c     # Horizontal velocity
      call annulus_covariant_basis(x, y, t1,t2) 

c     # Set non-zeros derivs only
      if (example .eq. 0) then
c        # Rigid body rotation        
         u(1) = revs_per_s
         u(2) = 0
      else
          if (example .eq. 1) then
              vcart(1) = cart_speed
              vcart(2) = 0
              vcart(3) = 0
          elseif (example .eq. 2) then
              A = amplitude
              tfinal = 0.25
              vcart(1) = cart_speed
              vcart(2) = pi2*A*cos(freq*pi2*t/tfinal)/tfinal;
          elseif (example .eq. 3) then
              A = amplitude
              tfinal = 0.25
              vcart(1) = cart_speed*pi*sin(pi*t/tfinal)/2.d0
              vcart(2) = 0
         elseif (example .eq. 4) then
             r = beta + (1-beta)*y
             th = theta(1) + (theta(2)-theta(1))*x
             w = 0.5
             vcart(1) = -(1-w)*pi2*r*sin(pi2*th)
             vcart(2) = (1+w)*pi2*r*cos(pi2*th)
             nc = sqrt(vcart(1)**2 + vcart(2)**2)

             vcart(1) = -vcart(1)/nc
             vcart(2) = -vcart(2)/nc
         endif

          vcart(3) = 0

          t1_dot_vcart = annulus_dot(t1,vcart)
          t2_dot_vcart = annulus_dot(t2,vcart)

          u(1) = t1_dot_vcart/t1_norm2         
          u(2) = t2_dot_vcart/t2_norm2
      endif

      end


      subroutine annulus_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta_vec(2)
      common /annulus_comm/ beta, theta_vec

      double precision r1, r1x, r1xx, r1y, r1yy, r1xy
      double precision R,  Rx,  Rxx,  Ry,  Ryy,  Rxy
      double precision f(3), fx(3), fy(3), fxx(3), fyy(3), fxy(3)
      double precision g(3), gx(3), gy(3), gxx(3), gyy(3), gxy(3)
      double precision t1(3), t2(3), a11, a12, a22, a21
      double precision det, a11inv, a22inv, a12inv, a21inv
      double precision pi4
      double precision annulus_dot


      integer k, kk, i
      logical compute_covariant, compute_contravariant
      logical compute_derivatives, b(32)

      double precision theta, thetax, thetay

      if (flag > 7) then
          write(6,*) 'psi.f : flag > 7'
          stop
      endif

c     # flag = 0      NA
c     # flag = 1      Covariant basis only
c     # flag = 2      NA
c     # flag = 3      Covariant + contravariant basis
c     # flag = 4      Derivatives only
c     # flag = 5      Derivatives + covariant basis
c     # flag = 6      NA
c     # flag = 7      Covariant + contravariant + derivatives


      do i = 1,bit_size(flag)
          b(i) = btest(flag,i-1)          
      enddo

      compute_covariant     = b(1) .or. b(2)
      compute_contravariant = b(2)
      compute_derivatives   = b(3)

      pi4 = pi2*pi2


c      r = beta + (1-beta)*yc
c      xp = r*cos(2*pi*xc)
c      yp = r*sin(2*pi*xc)

      theta = pi2*(theta_vec(1) + (theta_vec(2)-theta_vec(1))*x)
      thetax = pi2*(theta_vec(2)-theta_vec(1))
      thetay = 0

      r = beta + (1-beta)*y
      rx = 0
      ry = 1-beta
      if (compute_covariant) then
          t(1,1) = -thetax*r*sin(theta)
          t(2,1) = thetax*r*cos(theta)
          t(3,1) = 0

          t(1,2) = ry*cos(theta)
          t(2,2) = ry*sin(theta)
          t(3,2) = 0

      endif

      R     = beta + (1-beta)*y
      Rx    = 0
      Rxx   = 0
      Ry    = 1-beta
      Ryy   = 0
      Rxy   = 0

      f(1)  = cos(theta);
      fx(1) = -thetax*sin(theta)
      fy(1) = 0

      f(2)  = sin(theta);
      fx(2) = thetax*cos(theta)
      fy(2) = 0

      f(3) = 0
      fx(3) = 0
      fy(3) = 0

      g(1)  = R
      g(2)  = R
      g(3) = 0

      gx(1) = Rx
      gx(2) = Rx
      gx(3) = 0

      gy(1) = Ry
      gy(2) = Ry
      gy(3) = 0

      if (compute_covariant) then
          do k = 1,3
              t(k,1) = gx(k)*f(k) + g(k)*fx(k);
              t(k,2) = gy(k)*f(k) + g(k)*fy(k);
          enddo
      endif

      if (compute_contravariant) then
          do k = 1,3
              t1(k) = t(k,1)
              t2(k) = t(k,2)
          enddo

c         # Compute grad psi(xi,eta) 
          a11 = annulus_dot(t1,t1)
          a22 = annulus_dot(t2,t2)
          a12 = annulus_dot(t1,t2)
          a21 = a12

c         # Determinant
          det = a11*a22 - a12*a21

c         # Contravariant vectors
          a11inv = a22/det
          a22inv = a11/det
          a12inv = -a12/det
          a21inv = -a21/det     
          do k = 1,3
              tinv(k,1) = a11inv*t(k,1) + a12inv*t(k,2)
              tinv(k,2) = a21inv*t(k,1) + a22inv*t(k,2)
          end do
      endif

      if (compute_derivatives) then
          fxx(1) = -thetax**2*cos(theta)
          fyy(1) = 0
          fxy(1) = 0

          fxx(2) = -thetax**2*sin(theta)
          fyy(2) = 0
          fxy(2) = 0

          fxx(3) = 0
          fxy(3) = 0
          fyy(3) = 0

          gxx(1) = Rxx
          gxx(2) = Rxx
          gxx(3) = 0

          gyy(1) = Ryy
          gyy(2) = Ryy
          gyy(3) = 0

          gxy(1) = Rxy
          gxy(2) = Rxy
          gxy(3) = 0

          do k = 1,3
c             # d(t1)/dx = d(g*fx + gx*f)/dx
              tderivs(k,1,1) = g(k)*fxx(k) + 2*fx(k)*gx(k) + gxx(k)*f(k)

c             # d(t1)/dy = d(g*fx + gx*f)/dy       
              tderivs(k,1,2) = g(k)*fxy(k) + gy(k)*fx(k) + 
     &                    gx(k)*fy(k) + gxy(k)*f(k)

c             # d(t2)/dx = d(g*fy + gy*f)/dx       
              tderivs(k,2,1) = g(k)*fxy(k) + gx(k)*fy(k) + 
     &                    gy(k)*fx(k) + gxy(k)*f(k)

c             # d(t2)/dy = d(g*fy + gy*f)/dy         
              tderivs(k,2,2) = g(k)*fyy(k) + 2*fy(k)*gy(k) + gyy(k)*f(k)
          enddo
      endif

      end


