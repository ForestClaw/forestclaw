c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=torus; 1=twisted torus)
c     # ------------------------------------------------------------
      double precision function torus_psi(x,y)
      implicit none

      double precision x, y

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision psi

      if (beta .ne. 0) then
          write(6,'(A,A)') 'psi (psi.f) : Streamfunction only works ',
     &         'for beta == 0'
          stop
      endif

      psi = -pi2*revs_per_s*alpha*(pi2*y + alpha*sin(pi2*y))

      torus_psi = psi

      end


      subroutine torus_psi_derivs(x,y,p,px,py)
      implicit none

      double precision x, y
      double precision p, px, py

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha,beta

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      if (beta .ne. 0) then
          write(6,'(A,A)') 'psi (psi.f) : Streamfunction only works ',
     &         'for beta == 0'
          stop
      endif

c     # Stream function for rigid body rotation.      
      p = -(pi2*revs_per_s)*alpha*(pi2*y + alpha*sin(pi2*y))
      px = 0
      py = -(pi2)**2*revs_per_s*alpha*(1 + alpha*cos(pi2*y))

      end

      subroutine torus_velocity_components(x,y,u)
      implicit none

      double precision x, y, u(2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision s

      if (example .eq. 0) then
          u(1) = revs_per_s
          u(2) = 0
      elseif (example .eq. 1) then
c          s = sqrt(2.d0)
c          u(1) = s*cos(8*pi*x)
c          u(2) = s*sin(8*pi*y)   
           u(1) = 0
           u(2) = 1
      endif


      end

      subroutine torus_velocity_derivs(x,y,u,uderivs)
      implicit none

      double precision x, y, u(2)
      double precision uderivs(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision s, pim

      if (example .eq. 0) then
          u(1) = revs_per_s
          u(2) = 0
          uderivs(1) = 0
          uderivs(2) = 0
          uderivs(3) = 0
          uderivs(4) = 0
      elseif (example .eq. 1) then
c         s = sqrt(2.d0)
c         pim = 8*pi
c         u(1) = s*cos(pim*x)
c         u(2) = s*sin(pim*y)   
cc        # uderivs = [u1x u1y; u2x u2y]          
c         uderivs(1) = -s*pim*sin(pim*x)
c         uderivs(2) = 0;
c         uderivs(3) = 0; 
c         uderivs(4) = s*pim*cos(pim*y)
          u(1) = 0
          u(2) = 1
          uderivs(1) = 0
          uderivs(2) = 0
          uderivs(3) = 0
          uderivs(4) = 0           
      endif


      end


      subroutine torus_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision r1, r1x, r1xx, r1y, r1yy, r1xy
      double precision R,  Rx,  Rxx,  Ry,  Ryy,  Rxy
      double precision f(3), fx(3), fy(3), fxx(3), fyy(3), fxy(3)
      double precision g(3), gx(3), gy(3), gxx(3), gyy(3), gxy(3)
      double precision t1(3), t2(3), a11, a12, a22, a21
      double precision det, a11inv, a22inv, a12inv, a21inv
      double precision pi4
      double precision torus_dot

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

      r1    = alpha*(1 + beta*sin(pi2*x))
      r1x   = pi2*alpha*beta*cos(pi2*x)
      r1xx  = -pi4*alpha*beta*sin(pi2*x)

      R     = 1 +   r1*cos(pi2*y)
      Rx    =      r1x*cos(pi2*y)
      Rxx   =     r1xx*cos(pi2*y)
      Ry    =  -pi2*r1*sin(pi2*y)
      Ryy   =  -pi4*r1*cos(pi2*y)
      Rxy   = -pi2*r1x*sin(pi2*y)

      f(1)  = cos(pi2*x);
      fx(1) = -pi2*sin(pi2*x)
      fy(1) = 0;

      f(2)  = sin(pi2*x);
      fx(2) = pi2*cos(pi2*x)
      fy(2) = 0;

      f(3)  = sin(pi2*y);
      fx(3) = 0;
      fy(3) = pi2*cos(pi2*y);

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
          do k = 1,3
              t1(k) = t(k,1)
              t2(k) = t(k,2)
          enddo

c         # Compute grad psi(xi,eta) 
          a11 = torus_dot(t1,t1)
          a22 = torus_dot(t2,t2)
          a12 = torus_dot(t1,t2)
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
          fyy(1) = 0;
          fxy(1) = 0;

          fxx(2) = -pi4*sin(pi2*x)
          fyy(2) = 0;
          fxy(2) = 0;

          fxx(3) = 0;
          fxy(3) = 0;
          fyy(3) = -pi4*sin(pi2*y)

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

      subroutine torus_transform_coordinates(a,b,x,y,mapping)
      implicit none

      double precision a,b, x,y
      integer mapping
      double precision l0(4), l1(4), l2(4)

      data l0 /1., 0., 0., 1./
      data l1 /1., 0., 1., 1./
      data l2 /1., -1., 1., 0./

c     # This compute (x,y) from a1*t1 + a2*t2, where
c     # t1, t2 are the columns of the matrix L.

      if (mapping .eq. 0) then
          x = a
          y = b
      elseif (mapping .eq. 1)  then
          x = l1(1)*a + l1(2)*b
          y = l1(3)*a + l1(4)*b
      elseif (mapping .eq. 2) then
          x = l2(1)*a + l2(2)*b
          y = l2(3)*a + l2(4)*b
      endif

      end





