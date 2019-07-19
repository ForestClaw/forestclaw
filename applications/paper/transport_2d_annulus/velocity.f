c     # ------------------------------------------------------------
c     # Compute edge centered and cell-centered velocity fields
c     # 
c     # Edge : Compute velocity from a stream function.  
c     # This defines the average velocity at each edge
c     # 
c     #      u = curl \Psi  (div u = 0)
c     # 
c     # Center : Defines cell-centered velocities from a streamfunction
c     # 
c     #      u = n cross grad \Psi  (div u = 0)
c     # 
c     # or using basis functions
c     # 
c     #      u = u1*t1 + u2*t2   (div u might not be zero)
c     # 
c     # NOTE: All arguments to routines here should be in computational 
c     # coordinates (x,y) in [0,1]x[0,1].  Mapping from brick domains
c     # to [0,1]x[0,1] should be done from calling routines, as should 
c     # any mappings that convert from an orthogonal coordinate system
c     # to a non-orthogonal system.
c     # ------------------------------------------------------------


c     # ------------------------------------------------------------
c     # Center : u = n cross grad \Psi   (div u = 0)
c     # 
c     # Center : u = u1*tau1 + u2*tau2   (div u might not be zero)
c     # ------------------------------------------------------------
      subroutine annulus_center_velocity(x,y,t,vel)
      implicit none

      double precision x,y,vel(3), t

      double precision t1(3), t2(3), u(2)

      integer k

c     # Vector field defined as u1*tau1 + u2*tau2    

      call annulus_covariant_basis(x, y, t1,t2)
      call annulus_velocity_components(x,y,t,u)

      do k = 1,3
          vel(k) = u(1)*t1(k) + u(2)*t2(k)
      enddo

      end

      subroutine annulus_covariant_basis(x,y,t1,t2)
      implicit none

      double precision x,y,t1(3),t2(3)
      double precision t(3,2),tinv(3,2), tderivs(3,2,2)
      integer flag, k

c     # Compute covariant derivatives only
      flag = 1
      call annulus_basis_complete(x,y, t, tinv, tderivs, flag)

      do k = 1,3
          t1(k) = t(k,1)
          t2(k) = t(k,2)
      enddo

      end


      subroutine annulus_contravariant_basis(x,y,t1inv,t2inv)
      implicit none

      double precision x,y
      double precision t1inv(3), t2inv(3)
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      integer k, flag


      flag = 3
      call annulus_basis_complete(x,y, t, tinv,tderivs, flag)
                                          

      do k = 1,3
          t1inv(k) = tinv(k,1)
          t2inv(k) = tinv(k,2)
      end do

      end

c     # ----------------------------------------------------------------
c     # Utility functions
c     # ----------------------------------------------------------------

      subroutine annulus_cross(u,v,uxv,w)
      implicit none

      double precision u(3),v(3),uxv(3),w, annulus_dot
      integer k

      uxv(1) =   u(2)*v(3) - u(3)*v(2)
      uxv(2) = -(u(1)*v(3) - u(3)*v(1))
      uxv(3) =   u(1)*v(2) - u(2)*v(1)

      w = annulus_dot(uxv,uxv)
      w = sqrt(w)      

      end

      double precision function annulus_dot(u,v)
      implicit none

      double precision u(3),v(3)

      annulus_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      end





