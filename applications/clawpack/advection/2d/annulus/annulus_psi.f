c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=annulus; 1=twisted annulus)
c     # ------------------------------------------------------------
      double precision function annulus_psi(x,y)
      implicit none

      double precision x, y

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision r2

c     # annulus : Solid body rotation
      r2 = x**2 + y**2
      annulus_psi = -0.5d0*pi2*revs_per_s*r2

      end


      subroutine annulus_psi_derivs(x,y,p,px,py)
      implicit none

      double precision x, y
      double precision p, px, py

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision r, r2, rx, ry

c     # r2 = x**2 + y**2
c     # annulus_psi = 0.5d0*pi2*revs_per_s*r2

      r2 = x**2 + y**2
      r = sqrt(r)
      rx = x/r
      ry = y/r

      p = 0.5d0*pi2*revs_per_s*r2
      px = pi2*revs_per_s*r*rx
      py = pi2*revs_per_s*r*ry

      end

      subroutine annulus_velocity_components(x,y,u)
      implicit none

      double precision x, y, u(2)
      double precision uderivs(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision s

c     # uderivs not used here
      call annulus_velocity_derivs(x,y,u,uderivs)

      end

      subroutine annulus_velocity_derivs(x,y,u,uderivs)
      implicit none

      double precision x, y, u(2)
      double precision uderivs(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision s, pim, u1x, u1y, u2x, u2y
      integer k


c     # uderivs(1) = u1x      
c     # uderivs(2) = u1y      
c     # uderivs(3) = u2x      
c     # uderivs(4) = u2y      

      u1x = 0
      u1y = 0
      u2x = 0
      u2y = 0

c     # Set non-zeros derivs only
      s = sqrt(2.d0)
      if (example .eq. 0) then
c        # Rigid body rotation        
         u(1) = revs_per_s
         u(2) = 0
      elseif (example .eq. 1) then
c        # Conservative for all solvers (rp=1,2,3,4)               
         u(1) = s*(cos(pi*x)**2 + 0.5d0)         
         u(2) = s*(sin(pi*y)**2 + 0.5d0)
         u1x = -pi*s*sin(pi*x)
         u2y =  pi*s*cos(pi*y)
      elseif (example .eq. 2) then
         u(1) = s*(cos(pi*x)**2 - 0.5d0)
         u(2) = s*(sin(pi*y)**2 - 0.5d0)
         u1x = -pi*s*sin(pi*x)
         u2y =  pi*s*cos(pi*y)
      else
         write(6,'(A,A)') 'annulus_psi : ',
     &              'No valid example provided'
         stop
      endif

      uderivs(1) = u1x
      uderivs(2) = u1y
      uderivs(3) = u2x
      uderivs(4) = u2y

      end


      subroutine annulus_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta
      common /annulus_comm/ beta

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

      r = beta + (1-beta)*y
      rx = 0
      ry = 1-beta
      if (compute_covariant) then
          t(1,1) = -pi2*r*sin(pi2*x)
          t(2,1) = pi2*r*cos(pi2*x)
          t(3,1) = 0

          t(1,2) = ry*cos(pi2*x)
          t(2,2) = ry*sin(pi2*x)
          t(3,2) = 0

      endif

      R     = beta + (1-beta)*y
      Rx    = 0
      Rxx   = 0
      Ry    = 1-beta
      Ryy   = 0
      Rxy   = 0

      f(1)  = cos(pi2*x);
      fx(1) = -pi2*sin(pi2*x)
      fy(1) = 0

      f(2)  = sin(pi2*x);
      fx(2) = pi2*cos(pi2*x)
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
          fxx(1) = -pi4*cos(pi2*x)
          fyy(1) = 0
          fxy(1) = 0

          fxx(2) = -pi4*sin(pi2*x)
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

      subroutine annulus_transform_coordinates(a,b,x,y,mapping)
      implicit none

      double precision a,b, x,y
      integer mapping
      double precision l0(4), l1(4)

      data l0 /1., 0., 0., 1./
      data l1 /1., -0.2d0, 0., 1./

c     # This compute (x,y) from a1*t1 + a2*t2, where
c     # t1, t2 are the columns of the matrix L.

c      x = a
c      y = b

  
      if (mapping .eq. 0) then
          x = a
          y = b
      elseif (mapping .eq. 1) then
          x = a*l1(1) + b*l1(2)
          y = a*l1(3) + b*l1(4)
      endif

      end





