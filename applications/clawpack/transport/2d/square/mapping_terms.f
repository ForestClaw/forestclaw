c     # ------------------------------------------------------------
c     # 
c     # Mapping terms needed to compute exact solution
c     # 
c     # NOTE: All arguments to routines here should be in computational 
c     # coordinates (x,y) in [0,1]x[0,1].  Mapping from brick domains
c     # to [0,1]x[0,1] should be done from calling routines, as should 
c     # any mappings that convert from an orthogonal coordinate system
c     # to a non-orthogonal system.
c     # ------------------------------------------------------------


      subroutine map_covariant_basis(blockno,x,y,t1,t2)
      implicit none

      double precision x,y,t1(3),t2(3)
      integer blockno
      double precision t(3,2),tinv(3,2), tderivs(3,2,2)
      integer flag, k

      integer*8 cont, get_context

      cont = get_context()

c     # Compute covariant derivatives only
      flag = 1
      call fclaw2d_map_c2m_basis(cont, blockno, 
     &                 x,y, t, tinv, tderivs, flag)

      do k = 1,3
          t1(k) = t(k,1)
          t2(k) = t(k,2)
      enddo

      end


      subroutine map_contravariant_basis(blockno,x,y,t1inv,t2inv)
      implicit none

      double precision x,y
      integer blockno
      double precision t1inv(3), t2inv(3)
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      integer k, flag

      integer*8 cont, get_context

      cont = get_context()

      flag = 3
      call fclaw2d_map_c2m_basis(cont, blockno, 
     &             x,y, t, tinv,tderivs, flag)
                                          

      do k = 1,3
          t1inv(k) = tinv(k,1)
          t2inv(k) = tinv(k,2)
      end do

      end

      subroutine map_christoffel_sym(blockno,x,y,g) 
      implicit none

      double precision x, y
      integer blockno
      double precision g(2,2,2)

      double precision map_dot 

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision t(3,2), tinv(3,2), tderivs(3,2,2)
      double precision s(3,2), tij(3), sk(3)

      integer i,j,k, m, flag

      integer*8 cont, get_context

      cont = get_context()


c     # Compute covariant and derivatives
      flag = 7
      call fclaw2d_map_c2m_basis(cont, blockno, 
     &             x,y, t, tinv,tderivs, flag)
                                          

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

      double precision function map_divergence(blockno,x,y)
      implicit none

      double precision x,y
      integer blockno

      double precision u(2), uderivs(4), g(2,2,2)
      double precision D11, D22

c     # Get g(i,j,k), g = \Gamma(i,j,k)
      call velocity_derivs(blockno,x,y,u,uderivs)
      call map_christoffel_sym(blockno,x,y,g) 

      D11 = uderivs(1) + u(1)*g(1,1,1) + u(2)*g(1,2,1)
      D22 = uderivs(4) + u(1)*g(2,1,2) + u(2)*g(2,2,2)

      map_divergence = D11 + D22

      end


c     # ----------------------------------------------------------------
c     # Utility functions
c     # ----------------------------------------------------------------

      subroutine map_cross(u,v,uxv,w)
      implicit none

      double precision u(3),v(3),uxv(3),w, map_dot
      integer k

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





