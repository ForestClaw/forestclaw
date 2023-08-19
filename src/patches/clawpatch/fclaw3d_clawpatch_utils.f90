!> @file
!> Clawpatch utilities module

module fclaw3d_clawpatch_utils
implicit none
contains
      !     ---------------------------------------------------------------
      !>    @brief Determines the tranformation from a coarse patch to it's
      !>    half-sized neighbor's coordinate system
      !>
      !>    @param[in]  the pointer to the fclaw2d_patch_transform_data struct
      !>    @param[out] a the 2x2 tranform matrix for the i,j indexes of a patch
      !>    @param[out] f the transform vector
      !     ---------------------------------------------------------------
      subroutine build_transform(transform_ptr,a,f)
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
      logical function is_valid_average(i,j,k,mx,my,mz)
            implicit none

            integer :: i,j,k,mx,my,mz
            logical :: i1, j1, k1

            i1 = 1 .le. i .and. i .le. mx
            j1 = 1 .le. j .and. j .le. my
            k1 = 1 .le. k .and. k .le. mz

            is_valid_average = i1 .and. j1 .and. k1

      end function is_valid_average

      logical function is_valid_interp(i,j,k,mx,my,mz,mbc)
            implicit none

            integer i,j,k,mx,my,mz,mbc
            logical i1, j1, k1

            i1 = 1-mbc .le. i .and. i .le. mx+mbc
            j1 = 1-mbc .le. j .and. j .le. my+mbc
            k1 = 1-mbc .le. k .and. k .le. mz+mbc

            is_valid_interp = i1 .and. j1 .and. k1

      end function is_valid_interp

      subroutine interior_face_bounds( &
                  iface,mx,my,mz,mbc, &
                  i_start, j_start, k_start, &
                  i_end,   j_end,   k_end)
            implicit none
            integer mx, my, mz, mbc, iface
            integer i_start, j_start, k_start
            integer i_end,   j_end,   k_end
            integer axis
            logical upper

            axis = iface/2
            upper= btest(iface,0)
    
    
            if (axis .eq. 0) then
                if (upper) then
                    i_start = mx-mbc+1
                    i_end = mx
                else
                    i_start = 1
                    i_end = mbc
                endif
                j_start = 1
                j_end = my
                k_start = 1
                k_end = mz
            else if(axis .eq. 1) then
                i_start = 1
                i_end = mx
                if (upper) then
                    j_start = my-mbc+1
                    j_end = my
                else
                    j_start = 1
                    j_end = mbc
                endif
                k_start = 1
                k_end = mz
            else if(axis .eq. 2) then
                i_start = 1
                i_end = mx
                j_start = 1
                j_end = my
                if (upper) then
                    k_start = mz-mbc+1
                    k_end = mz
                else
                    k_start = 1
                    k_end = mbc
                endif
            endif
      end subroutine interior_face_bounds


end module fclaw3d_clawpatch_utils