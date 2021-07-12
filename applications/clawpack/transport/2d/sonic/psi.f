c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # Stream function depends on these global variables
c     #         example (0=incompressible; 1=compressible)
c     #         mapping (0=swirl; 1=twisted swirl)
c     # ------------------------------------------------------------
      double precision function swirl_psi(x,y)
      implicit none

      double precision x, y

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      integer example
      common /example_comm/ example  

      double precision psi, r

      r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

c     # Rigid body rotation
c     psi = r**2
 
c     # Filament formation (negative for clockwise rotation)
      psi = (4.d0/3.d0)*r**3

      swirl_psi = psi

      end


      subroutine swirl_psi_derivs(x,y,p,px,py)
      implicit none

      double precision x, y
      double precision p, px, py

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision r, rx, ry

      r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)    
      if (r .eq. 0) then
           write(5,*) 'swirl_psi_derivs (psi.f) : r .eq. 0'
           stop
      endif
      rx = x/r
      ry = y/r

      p = (4.d0/3.d0)*r**3          
      px = 4*r**2*rx
      py = 4*r**2*ry

      end

      subroutine swirl_velocity_components(x,y,u)
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

      call swirl_velocity_derivs(x,y,u,uderivs)

      end

      subroutine swirl_velocity_derivs(x,y,u,uderivs)
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
c        # Conservative for all solvers (rp=1,2,3,4)          
         u(1) = 1
         u(2) = 1
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
      else if (example .ge. 3) then         
c         u = 1
c         u = s*(sin(pi*yc)**2 + 0.5d0)
c         u = s*sin(pi*yc)

c         v = 0
c         v = -s*(sin(pi*xc)**2 + 0.5d0)
c         v = s*sin(pi*xc)

c         r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

c        # Rigid body rotation
c        psi = r**2
c         u(1) = 2*(yc-0.5)
c         u(2) = -2*(xc-0.5)
 
c        # Filament formation (negative for clockwise rotation)
c         psi = (4.d0/3.d0)*r**3

      else
         write(6,'(A,A)') 'clawpack46_setaux : ',
     &              'No valid example provided'
         stop
      endif

      uderivs(1) = u1x
      uderivs(2) = u1y
      uderivs(3) = u2x
      uderivs(4) = u2y

      end


      subroutine swirl_basis_complete(x,y, t, tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /swirl_comm/ alpha, beta

      double precision pi4
      double precision swirl_dot

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

      if (compute_covariant) then
          do k = 1,3
              t(k,1) = 0
              t(k,2) = 0
          enddo
          t(1,1) = 1
          t(2,2) = 1
      endif

      if (compute_contravariant) then
          do k = 1,3
              tinv(k,1) = 0
              tinv(k,2) = 0
          end do
          tinv(1,1) = 1
          tinv(2,2) = 1
      endif

      if (compute_derivatives) then
          do k = 1,3
              tderivs(k,1,1) = 0

c             # d(t1)/dy = d(g*fx + gx*f)/dy       
              tderivs(k,1,2) = 0

c             # d(t2)/dx = d(g*fy + gy*f)/dx       
              tderivs(k,2,1) = 0

c             # d(t2)/dy = d(g*fy + gy*f)/dy         
              tderivs(k,2,2) = 0
          enddo
      endif

      end

      subroutine swirl_transform_coordinates(a,b,x,y,mapping)
      implicit none

      double precision a,b, x,y
      integer mapping
      double precision l0(4), l1(4), l2(4)

      data l0 /1., 0., 0., 1./
      data l1 /1., 0., 1., 1./
      data l2 /1., -1., 1., 0./

c     # This compute (x,y) from a1*t1 + a2*t2, where
c     # t1, t2 are the columns of the matrix L.

      x = a
      y = b

  
c      if (mapping .eq. 0) then
c          x = a
c          y = b
c      elseif (mapping .eq. 1)  then
c          x = a*l1(1) + b*l1(2)
c          y = a*l1(3) + b*l1(4)
c      elseif (mapping .eq. 2) then
c          x = a*l2(1) + b*l2(2)
c          y = a*l2(3) + b*l2(4)
c      endif

      end





