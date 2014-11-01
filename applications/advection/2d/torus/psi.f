      double precision function psi(blockno,xc,yc,t)
      implicit none

      double precision xc, yc, t,r
      integer blockno
      double precision pi, r2, phi, pi2, alpha
      logical iscart, issphere
      double precision revs_per_s
      integer*8 cont, get_context

      double precision xp,yp,zp
      double precision xc1, yc1, zc1

      common /compi/ pi

      cont = get_context()

      call fclaw2d_map_brick2c(cont,
     &      blockno,xc,yc,xc1,yc1,zc1)

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      revs_per_s = 0.5d0

      alpha = 0.4d0
      pi2 = 2*pi
      if (iscart()) then
         psi = revs_per_s*(xc1 + pi*yc1/2.d0);
      elseif (issphere()) then
         psi = pi2*revs_per_s*zp
      else
c        # torus
         psi = revs_per_s*pi2*alpha*(pi2*yc1 + alpha*sin(pi2*yc1))
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
