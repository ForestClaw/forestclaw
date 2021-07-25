      subroutine torus_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      double precision r1, r1x, r1xx, r1y, r1yy, r1xy
      double precision R,  Rx,  Rxx,  Ry,  Ryy,  Rxy
      double precision f(3), fx(3), fy(3), fxx(3), fyy(3), fxy(3)
      double precision g(3), gx(3), gy(3), gxx(3), gyy(3), gxy(3)
c      double precision t1(3), t2(3), a11, a12, a22, a21
c      double precision det, a11inv, a22inv, a12inv, a21inv
      double precision pi4
      double precision map_dot

      double precision theta, phi
      double precision thetax, thetay, phix, phiy



      integer k, kk, i
      logical compute_covariant, compute_contravariant
      logical compute_derivatives, b(32)

      if (flag > 7) then
          write(6,*) 'torus_basis_complete : flag > 7'
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

c     # Torus mapping
c     #
c     #    r1 = alpha*(1 + beta*sin(pi2*x))
c     #    R = 1 + r1*cos(pi2*y)
c     #
c     #    xp =  R*cos(pi2*x)
c     #    yp =  R*sin(pi2*x)
c     #    zp = r1*sin(pi2*y)
c     #
c     #  Express X,Y,Z as 
c     # 
c     #     X(x,y) = g(x,y)*f(x,y) == g(1)*f(1)
c     #     Y(x,y) = g(x,y)*f(x,y) == g(2)*f(2)
c     #     Z(x,y) = g(x,y)*f(x,y) == g(3)*f(3)
c     # 
c     # Compute derivatives and use product rule to get higher order
c     # derivatives, i.e. 
c     #
c     #     dX/dx = g(1)*fx(1) + gx(1)*f(1)
c     #     dY/dx = g(2)*fx(2) + gx(2)*f(2)
c     #     dZ/dx = g(3)*fx(3) + gx(3)*f(3)
c     #

  
      call map_comp2toroidal_derivs(x,y,theta,phi,
     &           thetax, thetay, phix, phiy)


      r1    = alpha*(1 + beta*sin(theta))
      r1x   = alpha*beta*cos(theta)*thetax
      r1y   = 0
      r1xx  = -alpha*beta*sin(theta)*thetax*thetax
      r1yy  = 0
      r1xy  = 0

      R     = 1 +   r1*cos(phi)
      Rx    =      r1x*cos(phi) 
      Rxx   =     r1xx*cos(phi)
      Ry    =  -r1*sin(phi)*phiy
      Ryy   =  -r1*cos(phi)*phiy*phiy
      Rxy   = -r1x*sin(phi)*phiy

      f(1)  = cos(theta);
      fx(1) = -sin(theta)*thetax
      fy(1) = 0;

      f(2)  = sin(theta);
      fx(2) = cos(theta)*thetax
      fy(2) = 0;

      f(3)  = sin(phi);
      fx(3) = 0;
      fy(3) = cos(phi)*phiy

      g(1)  = R
      g(2)  = R
      g(3)  = r1

      gx(1) = Rx
      gx(2) = Rx
      gx(3) = r1x

      gy(1) = Ry
      gy(2) = Ry
      gy(3) = r1y

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
          fxx(1) = -cos(theta)*thetax**2
          fyy(1) = 0;
          fxy(1) = 0;

          fxx(2) = -sin(theta)*thetax**2
          fyy(2) = 0;
          fxy(2) = 0;

          fxx(3) = 0;
          fxy(3) = 0;
          fyy(3) = -sin(phi)*phiy**2

          gxx(1) = Rxx
          gxx(2) = Rxx
          gxx(3) = r1xx

          gyy(1) = Ryy
          gyy(2) = Ryy
          gyy(3) = r1yy

          gxy(1) = Rxy
          gxy(2) = Rxy
          gxy(3) = r1xy

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


      subroutine map_comp2torus(xc,yc,theta,phi)
      implicit none

      double precision xc,yc,theta, phi

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      double precision thetax, thetay, phix, phiy

c     # Map xc in [0,1] to theta in [-pi,pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      call map_comp2toroidal_derivs(xc,yc,theta,phi,
     &               thetax, thetay, phix, phiy)

      end

      subroutine map_comp2toroidal_derivs(xc,yc,theta,phi,
     &      thetax, thetay, phix, phiy)
      implicit none

      double precision xc,yc,theta, phi
      double precision  thetax, thetay, phix, phiy

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      double precision tr1, tr2, pr1, pr2


c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      


      tr1 = pi2*theta_range(1)
      tr2 = pi2*theta_range(2)

      theta = tr1 + (tr2-tr1)*xc
      thetax = tr2-tr1
      thetay = 0


      pr1 = pi2*phi_range(1)
      pr2 = pi2*phi_range(2)

      phi = pr1 + (pr2-pr1)*yc
      phix = 0
      phiy = pr2-pr1

      end

      subroutine map_toroidal2comp(theta,phi,xc,yc)
      implicit none

      double precision xc,yc,phi,theta

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      double precision tr1, tr2, pr1, pr2

c     # Map xc in [0,1] to theta in [0,2*pi]
c     # Map yc in [0,1] to phi in [-pi/2,pi/2]      

      tr1 = pi2*theta_range(1)
      tr2 = pi2*theta_range(2)

      pr1 = pi2*phi_range(1)
      pr2 = pi2*phi_range(2)

      xc = (theta-tr1)/(tr2-tr1)
      yc = (phi-pr1)/(pr2-pr1)

      end

      subroutine map2toroidal(xp,yp,zp,theta,phi)
      implicit none

      double precision xp,yp,zp,phi,theta

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta, theta_range(2), phi_range(2)
      common /torus_comm/ alpha, beta, theta_range, phi_range

      double precision r1, r

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

      r1 = alpha*(1 + beta*sin(theta))
      if (r1 .eq. 0) then
          write(6,*) 'torus_basis.f : r1 .eq. 0'
          stop
      endif

c     # Phi in [-pi/2, pi/2];  We need to transform it to 
c     # [0,2*pi]      
      phi = asin(zp/r1)
      r = sqrt(xp**2 + yp**2)
      if (r .le. 1) then
          if (phi .gt. 0) then
              phi = pi - phi
          else
              phi = -pi - phi
          endif
      endif
      if (phi < 0) then
          phi = phi + pi2
      endif



      end

      subroutine map2comp(xp,yp,zp,xc,yc)
      implicit none

      double precision xp,yp,zp, xc,yc

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision phi, theta

      call map2toroidal(xp,yp,zp,theta,phi)

      call map_toroidal2comp(theta,phi,xc,yc)

      end




cc     # Map to computational coordinates in [0,1]x[0,1]
c      subroutine map2comp(blockno,xc,yc,xp,yp,zp,xc1,yc1)
c      implicit none
c
c      double precision xc,yc,xp,yp,zp,xc1,yc1, zc1
c      integer blockno
c
c      integer*8 cont, get_context
c
c      cont = get_context()
c
cc     # We shouldn't be returning zc1 here - should always be zero in 2d
c      call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)
c
c
c      end

