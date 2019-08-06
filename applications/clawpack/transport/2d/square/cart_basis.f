      subroutine cart_basis_complete(blockno, x,y, t, 
     &                   tinv,tderivs, flag)
      implicit none
      
      double precision x,y 
      integer flag, blockno
      double precision t(3,2), tinv(3,2), tderivs(3,2,2)

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

c     # Needed to get Christoffel symbols
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

