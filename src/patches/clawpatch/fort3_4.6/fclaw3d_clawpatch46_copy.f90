!! # ----------------------------------------------------------
!! # Exchange routines - (i,j,mq) ordering
!! # ----------------------------------------------------------
!! # exchange_face_ghost
!! # exchange_corner_ghost
!! # exchange_phys_corner_ghost
!! # ----------------------------------------------------------


!! # Exchange edge ghost data with neighboring grid at same level.

subroutine fclaw3d_clawpatch46_fort_copy_face(mx,my,mz,mbc, & 
           meqn,qthis, qneighbor, iface, transform_ptr)

    implicit none

    integer mx,my,mz,mbc,meqn,iface
    integer*8 transform_ptr
    double precision     qthis(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer mq, k
    integer i,j

    integer a(2,2), f(2)

    !!call fclaw3d_clawpatch_build_transform_same(transform_ptr, a, f)

    if (iface .eq. 0) then
       do mq = 1,meqn
          do k = 1,mz
            do j = 1,my
               do i = 1-mbc,0
                  ! Lower side
                  qthis(i,j,k,mq) = qneighbor(i+mx,j,k,mq);
               enddo
            enddo
          enddo
       enddo
    else if (iface .eq. 1) then
       do mq = 1,meqn
          do k = 1,mz
            do j = 1,my
               do i = mx+1,mx+mbc
                  ! Upper side
                  qthis(i,j,k,mq) = qneighbor(i-mx,j,k,mq);
               enddo
            enddo
          enddo
       enddo
    else if (iface .eq. 2) then
       do mq = 1,meqn
          do k = 1,mz
            do j = 1-mbc,0
               do i = 1,mx
                  ! left side
                  qthis(i,j,k,mq) = qneighbor(i,j+my,k,mq);
               enddo
            enddo
          enddo
       enddo
    else if (iface .eq. 3) then
       do mq = 1,meqn
          do k = 1,mz
            do j = my+1,my+mbc
               do i = 1,mx
                  ! right side
                  qthis(i,j,k,mq) = qneighbor(i,j-my,k,mq);
               enddo
            enddo
          enddo
       enddo
    else if (iface .eq. 4) then
       do mq = 1,meqn
          do k = 1-mbc,0
            do j = 1,my
               do i = 1,mx
                  ! bottom
                  qthis(i,j,k,mq) = qneighbor(i,j,k+mz,mq);
               enddo
            enddo
          enddo
       enddo
    else if (iface .eq. 5) then
       do mq = 1,meqn
          do k = mz+1,mz+mbc
            do j = 1,my
               do i = 1,mx
                  ! top
                  qthis(i,j,k,mq) = qneighbor(i,j,k-mz,mq);
               enddo
            enddo
          enddo
       enddo
    endif



end subroutine fclaw3d_clawpatch46_fort_copy_face

subroutine fclaw3d_clawpatch46_fort_get_corner_start( &
            icorner,mx,my,mz,mbc, &
            i1, j1, k1)
    implicit none
    integer mx, my, mz, mbc, meqn, icorner
    integer i1, j1, k1

    if (icorner .eq. 0) then
         i1 = 1-mbc
         j1 = 1-mbc
         k1 = 1-mbc
    elseif (icorner .eq. 1) then
         i1 = mx+1
         j1 = 1-mbc
         k1 = 1-mbc
    elseif (icorner .eq. 2) then
         i1 = 1-mbc
         j1 = my+1
         k1 = 1-mbc
    elseif (icorner .eq. 3) then
         i1 = mx+1
         j1 = my+1
         k1 = 1-mbc
    elseif (icorner .eq. 4) then
         i1 = 1-mbc
         j1 = 1-mbc
         k1 = mz+1
    elseif (icorner .eq. 5) then
         i1 = mx+1
         j1 = 1-mbc
         k1 = mz+1
    elseif (icorner .eq. 6) then
         i1 = 1-mbc
         j1 = my+1
         k1 = mz+1
    elseif (icorner .eq. 7) then
         i1 = mx+1
         j1 = my+1
         k1 = mz+1
    endif
end subroutine fclaw3d_clawpatch46_fort_get_corner_start

subroutine fclaw3d_clawpatch46_fort_copy_corner(mx,my,mz,mbc,meqn, &
    qthis, qneighbor, this_icorner,transform_ptr)
    implicit none

    integer mx, my, mz, mbc, meqn, this_icorner
    integer*8 transform_ptr
    double precision     qthis(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer mq
    integer i_start, j_start, k_start
    integer i1,i2, j1,j2, k1,k2
    call fclaw3d_clawpatch46_fort_get_corner_start(this_icorner, &
                mx,my,mz,mbc, &
                i_start,j_start,k_start)

    !! # Do exchanges for all corners
    mq_loop : do mq = 1,meqn
         do k1 = k_start,k_start+mbc-1
             do j1 = j_start,j_start+mbc-1
                do i1 = i_start,i_start+mbc-1
                    call fclaw3d_clawpatch_transform_corner(i1,j1,k1,i2,j2,k2, transform_ptr)
                    qthis(i1,j1,k1,mq) = qneighbor(i2,j2,k2,mq)
                end do
            end do
        end do
    end do mq_loop

end subroutine fclaw3d_clawpatch46_fort_copy_corner
