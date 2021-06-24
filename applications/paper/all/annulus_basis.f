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

c      pi4 = pi2*pi2

      call map_comp2annulus_derivs(x,y,theta,r,
     &              thetax, thetay, rx, ry)      


c      r = beta + (1-beta)*yc
c      xp = r*cos(2*pi*xc)
c      yp = r*sin(2*pi*xc)



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
           call map_contravariant(t,tinv)          
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



      subroutine map_comp2annulus(xc,yc,theta,r)
      implicit none

      double precision xc,yc,theta, r

      double precision thetax, thetay, rx, ry

c     # Map xc in [0,1] to theta in [-pi,pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      call map_comp2annulus_derivs(xc,yc,theta,r,
     &               thetax, thetay, rx, ry)

      end

      subroutine map_comp2annulus_derivs(xc,yc,theta,r,
     &      thetax, thetay, rx, ry)
      implicit none

      double precision xc,yc,theta, r

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta_range(2)
      common /annulus_comm/ beta, theta_range

      double precision  thetax, thetay, rx, ry
      double precision tr1, tr2, pr1, pr2


c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      
  
      tr1 = pi2*theta_range(1)
      tr2 = pi2*theta_range(2)

      theta = tr1 + (tr2-tr1)*xc
      thetax = tr2-tr1
      thetay = 0

      r = beta + (1-beta)*yc
      rx = 0
      ry = 1-beta

      end

      subroutine map_annulus2comp(theta,r, xc,yc)
      implicit none

      double precision xc,yc,phi,theta, r

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta_range(2)
      common /annulus_comm/ beta, theta_range

      double precision tr1, tr2, pr1, pr2

c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      tr1 = pi2*theta_range(1)
      tr2 = pi2*theta_range(2)

      xc = (theta-tr1)/(tr2-tr1)
      yc = (r - beta)/(1-beta)

      end

      subroutine map2annulus(xp,yp,zp,theta,r)
      implicit none

      double precision xp,yp,zp,theta, r

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta_range(2)
      common /annulus_comm/ beta, theta_range

c      r1 = alpha*(1 + beta*sin(pi2*xc))
c      R = 1 + r1*cos(pi2*yc)
c
c      xp = R*cos(pi2*xc)
c      yp = R*sin(pi2*xc)
c      zp = r1*sin(pi2*yc)

  
      theta = atan2(yp,xp)    !! returns value in [-pi, pi]
      if (theta < 0) then
          theta = theta + pi2
      endif

      r = sqrt(xp**2 + yp**2)

      end

      subroutine map2comp(xp,yp,zp,xc,yc)
      implicit none

      double precision xp,yp,zp, xc,yc

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision r, theta

      call map2annulus(xp,yp,zp,theta,r)

      call map_annulus2comp(theta,r,xc,yc)

      end

