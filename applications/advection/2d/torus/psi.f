      double precision function psi(x,y,z,t)
      implicit none

      double precision x, y, z, t,r
      double precision pi, r2, phi, pi2, alpha
      logical iscart
      double precision revs_per_s

      common /compi/ pi

      revs_per_s = 1.d0

      alpha = 0.4d0
      pi2 = 2*pi
      if (iscart()) then
         psi = x
      else
c        # Finally works... !
         psi = revs_per_s*pi2*alpha*(pi2*y + alpha*sin(pi2*y))
      endif

      end


      subroutine get_vel(x,y,vvec,t)
      implicit none

      double precision x,y,vvec(3), t
      double precision cx,sx,cy,sy,r,pi
      double precision alpha, v, pi2
      double precision u_cart(3)
      double precision revs_per_s
      integer m

      common /compi/ pi

      alpha = 0.4d0
      pi2 = 2.d0*pi

c     # revs_per_s : revolutions/second
      revs_per_s = 1d0

      cx = cos(pi2*x)
      sx = sin(pi2*x)
      cy = cos(pi2*y)
      sy = sin(pi2*y)

      r = 1 + alpha*cy

c     # The Cartesian components of the velocity
      vvec(1) = -revs_per_s*r*pi2*sx
      vvec(2) = revs_per_s*r*pi2*cx
      vvec(3) = 0

      end

      double precision function torus_dot(u,v)
      implicit none
      double precision u(3),v(3)

      torus_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      end
