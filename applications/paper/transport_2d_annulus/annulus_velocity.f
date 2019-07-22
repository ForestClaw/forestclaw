c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=annulus; 1=twisted annulus)
c     # ------------------------------------------------------------



      subroutine annulus_center_velocity(x,y,t,vel)
      implicit none

      double precision x, y, vel(3), t

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s, vcart(2), amp, freq,cart_speed
      common /stream_comm/ revs_per_s, vcart, amp, freq, cart_speed

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer example
      common /example_comm/ example  

      double precision t1(3), t2(3), u(2)
c      double precision t1_norm2, t2_norm2
      double precision vc(3), annulus_dot
c      double precision t1_dot_vcart, t2_dot_vcart
      double precision xp,yp,zp, ravg, xc, d, tfinal, A
      double precision r, th, w, nc, vcr(2)

      integer k

      A = amp
      tfinal = 0.25


c     # Set non-zeros derivs only
      if (example .eq. 0) then
c         # Rigid body rotation        
          u(1) = revs_per_s
          u(2) = 0
          call annulus_covariant_basis(x, y, t1,t2) 
          do k = 1,3
            vel(k) = u(1)*t1(k) + u(2)*t2(k)
          end do
      else
          if (example .eq. 1) then
              vc(1) = vcart(1)
              vc(2) = vcart(2)
              vc(3) = 0
          elseif (example .eq. 2) then
              A = amp
              tfinal = 0.25
              vc(1) = cart_speed
              vc(2) = pi2*A*cos(freq*pi2*t/tfinal)/tfinal;
          elseif (example .eq. 3) then
              A = amp
              tfinal = 0.25
              vc(1) = cart_speed*pi*sin(pi*t/tfinal)/2.d0
              vc(2) = 0
         elseif (example .eq. 4) then
              A = amp
              tfinal = 0.25
              vc(1) = pi2*A*cos(freq*pi2*t/tfinal)/tfinal
c              vc(2) = cart_speed*pi*sin(pi*t/tfinal)/2.d0
              vc(2) = cart_speed
         elseif (example .eq. 5) then
             r = beta + (1-beta)*y
             th = theta(1) + (theta(2)-theta(1))*x
             w = amp
             vc(1) = -(1-w)*pi2*r*sin(pi2*th)
             vc(2) = (1+w)*pi2*r*cos(pi2*th)
             nc = sqrt(vc(1)**2 + vc(2)**2)

             vc(1) = -vc(1)/nc
             vc(2) = -vc(2)/nc
         endif
         vc(3) = 0

         do k = 1,3
            vel(k) = vc(k)
        end do

      endif
      end

      subroutine annulus_covariant_basis(x,y,t1,t2)
      implicit none

      double precision x,y,t1(3),t2(3)
      double precision theta, thetax, r, ry

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta_vec(2)
      common /annulus_comm/ beta, theta_vec

c      r = beta + (1-beta)*yc
c      xp = r*cos(2*pi*xc)
c      yp = r*sin(2*pi*xc)

      theta = pi2*(theta_vec(1) + (theta_vec(2)-theta_vec(1))*x)
      thetax = pi2*(theta_vec(2)-theta_vec(1))

      r = beta + (1-beta)*y
      ry = 1-beta

      t1(1) = -thetax*r*sin(theta)
      t1(2) = thetax*r*cos(theta)
      t1(3) = 0

      t2(1) = ry*cos(theta)
      t2(2) = ry*sin(theta)
      t2(3) = 0

      end

      double precision function annulus_dot(u,v)
      implicit none

      double precision u(3),v(3)

      annulus_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      end



