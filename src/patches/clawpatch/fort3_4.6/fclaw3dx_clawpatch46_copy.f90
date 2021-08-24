!! # ----------------------------------------------------------
!! # Exchange routines - (i,j,mq) ordering
!! # ----------------------------------------------------------
!! # exchange_face_ghost
!! # exchange_corner_ghost
!! # exchange_phys_corner_ghost
!! # ----------------------------------------------------------


!! # Exchange edge ghost data with neighboring grid at same level.

subroutine fclaw3dx_clawpatch46_fort_copy_face(mx,my,mz,mbc, & 
           meqn,qthis, qneighbor, iface, transform_ptr)

    implicit none

    integer mx,my,mz,mbc,meqn,iface
    integer*8 transform_ptr
    double precision     qthis(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer i,j,ibc,jbc,mq, idir, k
    integer i1,j1, i2, j2

    idir = iface/2

    !! # High side of 'qthis' exchanges with low side of
    !! # 'qneighbor'
    mq_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            if (idir .eq. 0) then
                do j = 1,my
                    do ibc = 1,mbc
                        !! # Exchange at low side of 'this' grid in
                        !! # x-direction (idir == 0)
                        if (iface .eq. 0) then
                            i1 = 1-ibc
                            j1 = j
                        elseif (iface .eq. 1) then
                            i1 = mx+ibc
                            j1 = j
                        endif
                        call fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,transform_ptr)
                        qthis(i1,j1,k,mq) = qneighbor(i2,j2,k,mq)
                    enddo
                enddo
            else
                do jbc = 1,mbc
                    do i = 1,mx
                        !! # Exchange at high side of 'this' grid in
                        !! # y-direction (idir == 1)
                        if (iface .eq. 2) then
                           i1 = i
                           j1 = 1-jbc
                        elseif (iface .eq. 3) then
                           i1 = i
                           j1 = my+jbc
                        endif
                        call fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,transform_ptr)
                        qthis(i1,j1,k,mq) = qneighbor(i2,j2,k,mq)
                    enddo
                enddo
            endif
        end do k_loop
    enddo mq_loop

end subroutine fclaw3dx_clawpatch46_fort_copy_face

subroutine fclaw3dx_clawpatch46_fort_copy_corner(mx,my,mz,mbc,meqn, &
    qthis, qneighbor, this_icorner,transform_ptr)
    implicit none

    integer mx, my, mz, mbc, meqn, this_icorner
    integer*8 transform_ptr
    double precision     qthis(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer mq, ibc, jbc
    integer i1, j1, i2, j2, k

    !! # Do exchanges for all corners
    mq_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            do ibc = 1,mbc
                do jbc = 1,mbc
                    if (this_icorner .eq. 0) then
                        i1 = 1-ibc
                        j1 = 1-jbc
                    elseif (this_icorner .eq. 1) then
                        i1 = mx+ibc
                        j1 = 1-jbc
                    elseif (this_icorner .eq. 2) then
                        i1 = 1 -ibc
                        j1 = my+jbc
                    else
                        i1 = mx+ibc
                        j1 = my+jbc
                    endif

                    !! # this routine is not yet complete, but the complete one
                    !! # can now be dropped in.
                    call fclaw2d_clawpatch_transform_corner(i1,j1,i2,j2, transform_ptr)
                    qthis(i1,j1,k,mq) = qneighbor(i2,j2,k,mq)
                end do
            end do
        end do k_loop
    end do mq_loop

end subroutine fclaw3dx_clawpatch46_fort_copy_corner
