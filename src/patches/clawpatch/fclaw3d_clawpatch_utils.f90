!> @file
!> Clawpatch utilities

!     ---------------------------------------------------------------
!>    @brief Determines the tranformation from a coarse patch to it's
!>    half-sized neighbor's coordinate system
!>
!>    @param[in]  the pointer to the fclaw2d_patch_transform_data struct
!>    @param[out] a the 2x2 tranform matrix for the i,j indexes of a patch
!>    @param[out] f the transform vector
!     ---------------------------------------------------------------
subroutine fclaw3d_clawpatch_build_transform(transform_ptr,a,f)
      implicit none

      integer :: a(3,3)
      integer*8 :: transform_ptr
      integer :: f(3)
      integer :: mi(8),mj(8),mk(8)
      integer :: i1,j1,k1

!     # Assume index mapping fclaw2d_transform_face_half has the
!     # the form
!     #
!     #       T(ic,jc) = A*(ic,jc) + F = (iff,jff)
!     #
!     # where (ic,jc) is the coarse grid index, and (iff,jff)
!     # is the first fine grid index.
!     #
!     # We can recover coefficients A(2,2) with the following
!     # calls to T.

      i1 = 0
      j1 = 0
      k1 = 0
      call fclaw3d_clawpatch_transform_face_half(i1,j1,k1,mi,mj,mk,transform_ptr)
      f(1) = mi(1)
      f(2) = mj(1)
      f(3) = mk(1)

      i1 = 1
      j1 = 0
      k1 = 0
      call fclaw3d_clawpatch_transform_face_half(i1,j1,k1,mi,mj,mk,transform_ptr)
      a(1,1) = mi(1) - f(1)
      a(2,1) = mj(1) - f(2)
      a(3,1) = mk(1) - f(3)

      i1 = 0
      j1 = 1
      k1 = 0
      call fclaw3d_clawpatch_transform_face_half(i1,j1,k1,mi,mj,mk,transform_ptr)
      a(1,2) = mi(1) - f(1)
      a(2,2) = mj(1) - f(2)
      a(3,2) = mk(1) - f(3)

      i1 = 0
      j1 = 0
      k1 = 1
      call fclaw3d_clawpatch_transform_face_half(i1,j1,k1,mi,mj,mk,transform_ptr)
      a(1,3) = mi(1) - f(1)
      a(2,3) = mj(1) - f(2)
      a(3,3) = mk(1) - f(3)

      end

      subroutine fclaw3d_clawpatch_build_transform_same(transform_ptr,a,f)
      implicit none

      integer a(2,2)
      integer*8 transform_ptr
      integer f(2)
      integer mi(4),mj(4)
      integer i1,j1

      i1 = 0
      j1 = 0
      call fclaw3d_clawpatch_transform_face(i1,j1,mi,mj,transform_ptr)
      f(1) = mi(1)
      f(2) = mj(1)

      i1 = 1
      j1 = 0
      call fclaw3d_clawpatch_transform_face(i1,j1,mi,mj,transform_ptr)
      a(1,1) = mi(1) - f(1)
      a(2,1) = mj(1) - f(2)

      i1 = 0
      j1 = 1
      call fclaw3d_clawpatch_transform_face(i1,j1,mi,mj,transform_ptr)
      a(1,2) = mi(1) - f(1)
      a(2,2) = mj(1) - f(2)

end
! --------------------------------------------------------------------
!> @brief checks if the index is valid for averaging
!>
!> @param[in] i, j, k the idnex to check
!> @param[in] mx, my, mz the number of cells in the x, y, z directions
!> @return true if the index is valid
! --------------------------------------------------------------------
logical function fclaw3d_clawpatch_is_valid_average(i,j,k,mx,my,mz)
      implicit none

      integer :: i,j,k,mx,my,mz
      logical :: i1, j1, k1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my
      k1 = 1 .le. k .and. k .le. mz

      fclaw3d_clawpatch_is_valid_average = i1 .and. j1 .and. k1

end function fclaw3d_clawpatch_is_valid_average

logical function fclaw3d_clawpatch_is_valid_interp(i,j,k,mx,my,mz,mbc)
      implicit none

      integer i,j,k,mx,my,mz,mbc
      logical i1, j1, k1

      i1 = 1-mbc .le. i .and. i .le. mx+mbc
      j1 = 1-mbc .le. j .and. j .le. my+mbc
      k1 = 1-mbc .le. k .and. k .le. mz+mbc

      fclaw3d_clawpatch_is_valid_interp = i1 .and. j1 .and. k1

end function fclaw3d_clawpatch_is_valid_interp