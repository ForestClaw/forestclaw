      subroutine sphere_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision f(3), fx(3), fy(3), fxx(3), fyy(3), fxy(3)
      double precision g(3), gx(3), gy(3), gxx(3), gyy(3), gxy(3)
      double precision map_dot


      integer k, kk, i
      logical compute_covariant, compute_contravariant
      logical compute_derivatives, b(32)

      double precision theta, thetax, thetay
      double precision phi, phix, phiy

      if (flag > 7) then
          write(6,*) 'sphere_basis.f : flag > 7'
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

c     # ---------------------------
c     # Spherical coordinates
c     # ---------------------------
c         xp = cos(phi)*cos(theta)
c         yp = cos(phi)*sin(theta)
c         zp = sin(phi)
c     # ---------------------------

 
cc     # Longitude
c      theta = pi2*x
c      thetax = pi2
c      thetay = 0
c
cc     # Latitude
c      phi = -pi/2.d0 + pi*y
c      phix = 0
c      phiy = pi

      call map_comp2angles_derivs(x,y,phi,theta,
     &                     thetax, thetay, phix, phiy)


      if (compute_covariant) then
c          t(1,1) = -thetax*cos(phi)*sin(theta)
c          t(2,1) = thetax*cos(phi)*cos(theta)
c          t(3,1) = 0

c          t(1,2) = -phiy*sin(phi)*cos(theta)
c          t(2,2) = -phiy*sin(phi)*sin(theta)
c          t(3,2) = phiy*cos(phi)          

      endif

c     # Express T(x,y) = g(y)*f(x)
c     #
c     #    T1(x,y) = g1(y)*f1(x)
c     #    T2(x,y) = g2(y)*f2(x)
c     #    T3(x,y) = g3(y)*f3(x)
c     #

      f(1)  = cos(theta)
      f(2)  = sin(theta)
      f(3) = 1

      fx(1) = -thetax*sin(theta)
      fx(2) = thetax*cos(theta)
      fx(3) = 0

      fy(1) = 0
      fy(2) = 0
      fy(3) = 0

      g(1)  = cos(phi)
      g(2)  = cos(phi)
      g(3) = sin(phi)

      gx(1) = 0
      gx(2) = 0
      gx(3) = 0

      gy(1) = -phiy*sin(phi)
      gy(2) = -phiy*sin(phi)
      gy(3) = phiy*cos(phi)

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
          fxx(2) = -thetax**2*sin(theta)
          fxx(3) = 0

          fyy(1) = 0
          fyy(2) = 0
          fyy(3) = 0

          fxy(1) = 0
          fxy(2) = 0
          fxy(3) = 0

          gxx(1) = 0
          gxx(2) = 0
          gxx(3) = 0

          gyy(1) = -phiy**2*cos(phi)
          gyy(2) = -phiy**2*cos(phi)
          gyy(3) = 0

          gxy(1) = 0
          gxy(2) = 0
          gxy(3) = 0

          do k = 1,3
c             # d(t1)/dx = d(g*fx + gx*f)/dx
              tderivs(k,1,1) = g(k)*fxx(k) + 
     &                     2*fx(k)*gx(k) + gxx(k)*f(k)

c             # d(t1)/dy = d(g*fx + gx*f)/dy       
              tderivs(k,1,2) = g(k)*fxy(k) + gy(k)*fx(k) + 
     &                    gx(k)*fy(k) + gxy(k)*f(k)

c             # d(t2)/dx = d(g*fy + gy*f)/dx       
              tderivs(k,2,1) = g(k)*fxy(k) + gx(k)*fy(k) + 
     &                    gy(k)*fx(k) + gxy(k)*f(k)

c             # d(t2)/dy = d(g*fy + gy*f)/dy         
              tderivs(k,2,2) = g(k)*fyy(k) + 
     &                    2*fy(k)*gy(k) + gyy(k)*f(k)
          enddo
      endif

      end



      subroutine map_comp2angles(xc,yc,phi,theta)
      implicit none

      double precision xc,yc,theta, thetax, thetay
      double precision phi, phix, phiy

      double precision pi, pi2
      common /compi/ pi, pi2

c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      call map_comp2angles_derivs(xc,yc,phi,theta,
     &               thetax, thetay, phix, phiy)

      end

      subroutine map_comp2angles_derivs(xc,yc,phi,theta,
     &      thetax, thetay, phix, phiy)
      implicit none

      double precision xc,yc,theta, thetax, thetay
      double precision phi, phix, phiy

      double precision pi, pi2
      common /compi/ pi, pi2

c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      theta = pi2*xc
      thetax = pi2
      thetay= 0

      phi = -pi/2 + pi*yc
      phix = 0
      phiy = pi

      end

      subroutine map_angles2comp(phi,theta,xc,yc)
      implicit none

      double precision xc,yc,phi,theta

      double precision pi, pi2
      common /compi/ pi, pi2


c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      xc = theta/pi2
      yc = (phi + pi/2)/pi

      end

      subroutine map2polar(xp,yp,zp,phi,theta)
      implicit none

      double precision xp,yp,zp,phi,theta

      double precision pi, pi2
      common /compi/ pi, pi2


      phi = pi/2 - acos(zp) 
      theta = atan2(yp,zp)
      if (theta < 0) then
          theta = theta + pi2
      endif


      end

      subroutine map2comp(xp,yp,zp,xc,yc)
      implicit none

      double precision xp,yp,zp, xc,yc

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision phi, theta

      call map2polar(xp,yp,zp,phi,theta)

      call map_angles2comp(phi,theta,xc,yc)

      end


