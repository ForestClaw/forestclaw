c     # ------------------------------------------------------------
c     # 
c     # Mapping terms needed to compute exact solution
c     # 
c     # These routines all rely on a basis function that returns 
c     # covariant and contravariant vectors for the domain.  For 
c     # example, the square domain uses basis vectors (1,0) and (0,1).
c     # whereas the annular domain uses the usual polar coordinate basis
c     # vectors.
c     # 
c     # NOTE : The blockno should only be used if the underlying domain
c     # cannot be described by a single Cartesian domain.  The 
c     # cubedsphere for example, is a multiblock description of 
c     # the spherical coordinate system.  But the spherical coordinate 
c     # system can be described with a single Cartesian block.
c     # ------------------------------------------------------------


      subroutine map_covariant_basis(x,y,t1,t2)
      implicit none

      double precision x,y,t1(3),t2(3)
      double precision t(3,2),tinv(3,2), tderivs(3,2,2)
      integer flag, k

      integer*8 cont, fclaw_map_get_context

      cont = fclaw_map_get_context()

c     # Compute covariant derivatives only
      flag = 1
      call fclaw_map_2d_c2m_basis(cont, 
     &                 x,y, t, tinv, tderivs, flag)

      do k = 1,3
          t1(k) = t(k,1)
          t2(k) = t(k,2)
      enddo

      end


      subroutine map_contravariant_basis(x,y,t1inv,t2inv)
      implicit none

      double precision x,y
      double precision t1inv(3), t2inv(3)
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      integer k, flag

      integer*8 cont, fclaw_map_get_context

      cont = fclaw_map_get_context()

      flag = 3
      call fclaw_map_2d_c2m_basis(cont, x,y, t, tinv,tderivs, flag)
                                          

      do k = 1,3
          t1inv(k) = tinv(k,1)
          t2inv(k) = tinv(k,2)
      end do

      end

      subroutine map_christoffel_sym(x,y,g) 
      implicit none

      double precision x, y
      double precision g(2,2,2)

      double precision map_dot 

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      double precision s(3,2), tij(3), sk(3)

      integer i,j,k, m, flag

      integer*8 cont, fclaw_map_get_context

      cont = fclaw_map_get_context()


c     # Compute covariant and derivatives
      flag = 7
      call fclaw_map_2d_c2m_basis(cont, x,y, t,tinv,tderivs,flag)
                                          

      do k = 1,3
          s(k,1) = tinv(k,1)
          s(k,2) = tinv(k,2)
      end do

      do i = 1,2
          do j = 1,2
              do k = 1,2
                  do m = 1,3
                      tij(m) = tderivs(m,i,j)
                      sk(m) = s(m,k) 
                  enddo
                  g(i,j,k) = map_dot(tij,sk)
              enddo
          enddo
      enddo

      end

      double precision function map_divergence(x,y, t)
      implicit none

      double precision x,y, t

      double precision u(2), vcart(3), derivs(4)
      double precision D11, D22, g(2,2,2)
      integer flag

c     # Get g(i,j,k), g = \Gamma(i,j,k)
      call velocity_derivs(x,y,t, u,vcart,derivs,flag)

      if (flag .eq. 0) then
c         # Velocity and derivatives are given in 
c         # spherical components       
          call map_christoffel_sym(x,y,g) 

          D11 = derivs(1) + u(1)*g(1,1,1) + u(2)*g(1,2,1)
          D22 = derivs(4) + u(1)*g(2,1,2) + u(2)*g(2,2,2)

          map_divergence = D11 + D22
      else
c         # Velocity and derivatives are given in 
c         # Cartesian components        
          map_divergence = derivs(1) + derivs(2) + derivs(3)
      endif

      end


c     # ----------------------------------------------------------------
c     # Utility functions
c     # ----------------------------------------------------------------

      subroutine map_cross(u,v,uxv,w)
      implicit none

      double precision u(3),v(3),uxv(3),w, map_dot

      uxv(1) =   u(2)*v(3) - u(3)*v(2)
      uxv(2) = -(u(1)*v(3) - u(3)*v(1))
      uxv(3) =   u(1)*v(2) - u(2)*v(1)

      w = map_dot(uxv,uxv)
      w = sqrt(w)      

      end

      double precision function map_dot(u,v)
      implicit none

      double precision u(3),v(3)

      map_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      end


c     # Compute contravariant basis vectors from the 
c     # covariant basis.
      subroutine map_contravariant(t,tinv)
      implicit none

      double precision t(3,2), tinv(3,2)

      double precision t1(3), t2(3)
      integer k

      double precision a11, a22, a12, a21, det, map_dot
      double precision a11inv, a22inv, a12inv, a21inv

      do k = 1,3
          t1(k) = t(k,1)
          t2(k) = t(k,2)
      enddo

c     # Compute grad psi(xi,eta) 
      a11 = map_dot(t1,t1)
      a22 = map_dot(t2,t2)
      a12 = map_dot(t1,t2)
      a21 = a12

c     # Determinant
      det = a11*a22 - a12*a21

c     # Contravariant vectors
      a11inv = a22/det
      a22inv = a11/det
      a12inv = -a12/det
      a21inv = -a21/det     
      do k = 1,3
        tinv(k,1) = a11inv*t(k,1) + a12inv*t(k,2)
        tinv(k,2) = a21inv*t(k,1) + a22inv*t(k,2)
      end do

      end

      subroutine map_diff_normalized(x,y,u, uderivs, 
     &         uderivs_norm)
      implicit  none

      double precision x,y,u(2), uderivs(4), uderivs_norm(4)

      integer*8 cont, fclaw_map_get_context

      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      double precision t1(3), t2(3), t1n2, t2n2
      double precision t1x(3), t1y(3), t2x(3), t2y(3)
      double precision map_dot, map_quotient_rule
      double precision t1n, t2n, dt1ndx, dt1ndy,dt2ndx,dt2ndy
      integer flag, k

      cont = fclaw_map_get_context()

      flag = 7
      call fclaw_map_2d_c2m_basis(cont,x,y,t,tinv, tderivs,flag)

      do k = 1,3
          t1(k) = t(k,1)
          t2(k) = t(k,2)
          t1x(k) = tderivs(k,1,1)
          t1y(k) = tderivs(k,1,2)
          t2x(k) = tderivs(k,2,1)
          t2y(k) = tderivs(k,2,2)
      enddo

      t1n2 = map_dot(t1,t1)
      t2n2 = map_dot(t2,t2)

      t1n = sqrt(t1n2)
      t2n = sqrt(t2n2)


c     # d(norm(t1))/dx = d(sqrt(t1 dot t1))/dx      
      dt1ndx = map_dot(t1,t1x)/t1n
      dt1ndy = map_dot(t1,t1y)/t1n

      dt2ndx = map_dot(t2,t2x)/t2n
      dt2ndy = map_dot(t2,t2y)/t2n


c     # d(u1/t1n)/dx      
c      uderivs_norm(1) = (t1n2*uderivs(1) - u(1)*map_dot(t1x,t1))/t1n3
      uderivs_norm(1) = map_quotient_rule(u(1),t1n,uderivs(1),dt1ndx)

c     # d(u1/t1n)/dy      
c      uderivs_norm(2) = (t1n2*uderivs(2) - u(1)*map_dot(t1y,t1))/t1n3
      uderivs_norm(2) = map_quotient_rule(u(1),t1n,uderivs(2),dt1ndy)

c     # d(u2/t2n)/dx      
c      uderivs_norm(3) = (t2n2*uderivs(3) - u(2)*map_dot(t2x,t2))/t2n3
      uderivs_norm(3) = map_quotient_rule(u(2),t2n,uderivs(3),dt2ndx)

c     # d(u2/t2n)/dy
c      uderivs_norm(4) = (t2n2*uderivs(4) - u(2)*map_dot(t2y,t2))/t2n3
      uderivs_norm(4) = map_quotient_rule(u(2),t2n,uderivs(4),dt2ndy)



      end


      double precision function map_quotient_rule(f,g,fx,gx)
      implicit none

      double precision f,g,fx,gx

c     # d(f/g)/dx = (g*fx - f*gx)/g^2
      map_quotient_rule = (g*fx - f*gx)/g**2

      end




