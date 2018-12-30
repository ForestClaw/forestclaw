      double precision function torus_psi(blockno,xc,yc,t)
      implicit none

      double precision xc, yc, t
      integer blockno
      double precision pi, alpha
      double precision revs_per_s
      integer*8 cont, get_context

      double precision xc1, yc1, zc1, pi2
      integer example, mapping

      common /compi/ pi
      common /torus_comm/ alpha, revs_per_s
      common /mapping_comm/ mapping

      double precision psi

      cont = get_context()

c     # This is not the torus mapping, but rather maps the brick to
c     # a unit square      
      call fclaw2d_map_brick2c(cont,
     &      blockno,xc,yc,xc1,yc1,zc1)

      pi2 = 2*pi

c     # Velocity is described in terms of (xi, eta)=(xc1,yc1) coordinates
      if (mapping .eq. 0) then
c        # Rigid body rotation
         psi = (pi2*revs_per_s)*alpha*
     &         (pi2*yc1 + alpha*sin(pi2*yc1))
      elseif (mapping .eq. 1) then
c        # Twisted torus stream function (to be used with usual torus map)
         psi = (pi2*revs_per_s)*alpha*
     &         (pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
      endif
      torus_psi = psi

      end

      subroutine torus_edge_velocity(blockno,xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(2),xd2(2), ds, vn, torus_psi,t
      double precision t1, t2
      integer blockno

      vn = (torus_psi(blockno,xd1(1),xd1(2),t) -
     &      torus_psi(blockno,xd2(1),xd2(2),t))/ds

      end

      subroutine torus_center_velocity(blockno,xc,yc,t, 
     &                                  nv, u, v, w)
      implicit none

      integer blockno
      double precision xc,yc,t,nv(3),kappa(3), u,v, w


      double precision pi, alpha
      double precision revs_per_s
      integer*8 cont, get_context

      double precision xc1, yc1, zc1, pi2
      double precision tau1(3), tau2(3)
      double precision R, Reta, Rxi, vel(3), uc, vc
      double precision a11, a22, a12, a21
      double precision a11inv, a22inv, a12inv, a21inv
      double precision tau1inv(3), tau2inv(3)
      double precision psi_xi, psi_eta, gradpsi(3), nvec(3)
      double precision vdotn, det, h1, h2
      integer k

      integer example, mapping
      common /compi/ pi
      common /torus_comm/ alpha, revs_per_s
      common /mapping_comm/ mapping

      cont = get_context()

c     # This is not the torus mapping, but rather maps the brick to
c     # a unit square      
      call fclaw2d_map_brick2c(cont,
     &      blockno,xc,yc,xc1,yc1,zc1)

      pi2 = 2*pi

      if (mapping .eq. 0) then
c        # Rigid body rotation (\Psi_y, -\Psi_x)
          psi_xi = 0
          psi_eta = (pi2)**2*revs_per_s*alpha*(1 + alpha*cos(pi2*yc1))

c         # Coordinate normals          
          R = 1 + alpha*cos(pi2*yc1)
          Reta = -pi2*alpha*sin(pi2*yc1)

c         # tau1
          tau1(1) = -pi2*R*sin(pi2*xc1)
          tau1(2) = pi2*R*cos(pi2*xc1)
          tau1(3) = 0

c         # tau2
          tau2(1) = Reta*cos(pi2*xc1)
          tau2(2) = Reta*sin(pi2*xc1)
          tau2(3) = pi2*alpha*cos(pi2*yc1)               

      elseif (mapping .eq. 1) then
c        # Twisted torus stream function (to be used with usual torus map)
c         psi = (pi2*revs_per_s)*alpha*(pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))

          psi_xi = (pi2)**2*revs_per_s*alpha*(1 + 
     &               alpha*cos(pi2*(xc1+yc1)))
          psi_eta = (pi2)**2*revs_per_s*alpha*(1 + 
     &               alpha*cos(pi2*(xc1+yc1)))

c         # Coordinate normals          
          R    = 1 +  alpha*cos(pi2*(xc1 + yc1))
          Rxi  = -pi2*alpha*sin(pi2*(xc1 + yc1))
          Reta = -pi2*alpha*sin(pi2*(xc1 + yc1))

          tau1(1) = Rxi*cos(pi2*xc1) - pi2*R*sin(pi2*xc1)
          tau1(2) = Rxi*sin(pi2*xc1) + pi2*R*cos(pi2*xc1)
          tau1(3) = pi2*alpha*cos(pi2*(xc1 + yc1))

          tau2(1) = Reta*cos(pi2*xc1)
          tau2(2) = Reta*sin(pi2*xc1)
          tau2(3) = pi2*alpha*cos(pi2*(xc1+yc1))
      endif


      a11 = 0
      a22 = 0
      a12 = 0
      do k = 1,3
          a11 = a11 + tau1(k)*tau1(k)
          a22 = a22 + tau2(k)*tau2(k)
          a12 = a12 + tau1(k)*tau2(k)
      end do

      det = a11*a22 - a12*a12 
      if (det .eq. 0) then
          write(6,*) 'velocity : Determinant is 0'
          stop
      endif

      a11inv = a22/det
      a22inv = a11/det
      a12inv = -a12/det
      a21inv = -a12/det     
      do k = 1,3
          tau1inv(k) = a11inv*tau1(k) + a12inv*tau2(k)
          tau2inv(k) = a21inv*tau1(k) + a22inv*tau2(k)
      end do

      do k = 1,3
          gradpsi(k) = psi_xi*tau1inv(k) + psi_eta*tau2inv(k)
      end do

      call torus_cross(tau1,tau2,nvec,w)

c     # Normalize surface normal
      do k = 1,3
          nvec(k) = nvec(k)/w
      end do

      call torus_cross(gradpsi,nvec,vel,w)

c     # Compute velocity in physical space uvec = u*tau1 + v*tau2
      vdotn = 0
      do k = 1,3
c         # Compute vdotn = v dot (surface normal)          
          vdotn = vdotn + vel(k)*nv(k)
      enddo


c     # Subtract out components of the surface normal
      do k = 1,3
          vel(k) = vel(k) - vdotn*nv(k)
      enddo

c     # Velocity field should satisfy (u,v,w) dot n = 0 so that we 
c     # get conservation.      
      u = vel(1)
      v = vel(2)
      w = vel(3)

      end

      subroutine torus_cross(u,v,uxv,w)
      implicit none

      double precision u(3),v(3),uxv(3),w
      integer k

      uxv(1) =   u(2)*v(3) - u(3)*v(2)
      uxv(2) = -(u(1)*v(3) - u(3)*v(1))
      uxv(3) =   u(1)*v(2) - u(2)*v(1)

      w = 0
      do k = 1,3
         w = w + uxv(k)*uxv(k)
      enddo
      w = sqrt(w)      

      end




