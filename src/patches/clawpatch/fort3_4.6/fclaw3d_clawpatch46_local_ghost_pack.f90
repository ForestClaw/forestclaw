subroutine fclaw3d_clawpatch46_fort_local_ghost_pack (mx,my,mz, mbc, & 
    meqn, mint,qdata,volume,qpack,psize, packmode,ierror)

    implicit none
    integer :: mx,my,mz, mbc,meqn,psize, mint
    integer :: packmode, ierror
    double precision :: qdata(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc + 1)
    double precision :: qpack(psize)

    integer :: packq, unpackq, packarea, unpackarea
    parameter(packq = 0, unpackq = 1, packarea = 2, unpackarea = 3)

    integer :: i,j,k,mq,ibc,jbc,kbc, count, count_final
    integer :: nghost
    logical :: packdata

    ierror = 0
    if (packmode .ne. packq .and. packmode .ne. unpackq .and. &
           packmode .ne. packarea .and. packmode .ne. unpackarea) then
        ierror = 1
        return
    endif

    packdata = packmode .eq. packq .or. packmode .eq. packarea

    nghost = mbc
    count = 1
    meqn_loop : do mq = 1,meqn
        !! # Face 0
        do k = 1-nghost,mz+nghost
            do j = 1-nghost,my+nghost
                do ibc = 1-nghost,mint
                    if (packdata) then
                        qpack(count) = qdata(ibc,j,k,mq)
                    else
                        qdata(ibc,j,k,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

        !! # Face 1
        do k = 1-nghost,mz+nghost
            do j = 1-nghost,my+nghost
                do ibc = mx-mint+1,mx+nghost
                    if (packdata) then
                        qpack(count) = qdata(ibc,j,k,mq)
                    else
                        qdata(ibc,j,k,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

        !! # Face 2
        do k = 1-nghost,mz+nghost
            do jbc = 1-nghost,mint
                do i = mint+1,mx-mint
                    if (packdata) then
                        qpack(count) = qdata(i,jbc,k,mq)
                    else
                        qdata(i,jbc,k,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

        !! # Face 3
        do k = 1-nghost,mz+nghost
            do jbc = my-mint+1,my+nghost
                do i = mint+1,mx-mint
                    if (packdata) then
                        qpack(count) = qdata(i,jbc,k,mq)
                    else
                        qdata(i,jbc,k,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

        !! # Face 4
        do kbc = 1-nghost,mint
            do j = mint+1,my-mint
                do i = mint+1,mx-mint
                    if (packdata) then
                        qpack(count) = qdata(i,j,kbc,mq)
                    else
                        qdata(i,j,kbc,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

        !! # Face 5
        do kbc = mz-mint+1,mz+nghost
            do j = mint+1,my-mint
                do i = mint+1,mx-mint
                    if (packdata) then
                        qpack(count) = qdata(i,j,kbc,mq)
                    else
                        qdata(i,j,kbc,mq) = qpack(count)
                    endif
                    count = count + 1
                end do
            end do
        end do

    end do meqn_loop

    count_final = count - 1
    if (packmode .ne. packarea .and. packmode .ne. unpackarea) then
        if (count_final .ne. psize) then
            write(6,*) 'Before volume packing/unpacking'
            write(6,*) 'count_final = ',count_final
            write(6,*) 'psize = ',psize
            ierror = 2
        endif
        return
    endif

    !! # Face 0
    do k = 1-nghost,mz+nghost
        do j = 1-nghost,my+nghost
            do ibc = 1-nghost,mint
                if (packdata) then
                    qpack(count) = volume(ibc,j,k)
                else
                    volume(ibc,j,k) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    !! # Face 1
    do k = 1-nghost,mz+nghost
        do j = 1-nghost,my+nghost
            do ibc = mx-mint+1,mx+nghost
                if (packdata) then
                    qpack(count) = volume(ibc,j,k)
                else
                    volume(ibc,j,k) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    !! # Face 2
    do k = 1-nghost,mz+nghost
        do jbc = 1-nghost,mint
            do i = mint+1,mx-mint
                if (packdata) then
                    qpack(count) = volume(i,jbc,k)
                else
                    volume(i,jbc,k) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    !! # Face 3
    do k = 1-nghost,mz+nghost
        do jbc = my-mint+1,my+nghost
            do i = mint+1,mx-mint
                if (packdata) then
                    qpack(count) = volume(i,jbc,k)
                else
                    volume(i,jbc,k) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    !! # Face 4
    do kbc = 1-nghost,mint
        do j = mint+1,my-mint
            do i = mint+1,mx-mint
                if (packdata) then
                    qpack(count) = volume(i,j,kbc)
                else
                    volume(i,j,kbc) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    !! # Face 5
    do kbc = mz-mint+1,mz+nghost
        do j = mint+1,my-mint
            do i = mint+1,mx-mint
                if (packdata) then
                    qpack(count) = volume(i,j,kbc)
                else
                    volume(i,j,kbc) = qpack(count)
                endif
                count = count + 1
            end do
        end do
    end do

    count_final = count-1
    if (count_final .ne. psize) then
        write(6,*) 'After volume packing/unpacking'
        write(6,*) 'psize = ', psize
        write(6,*) 'count_final = ', count_final
       ierror = 2
   endif

end subroutine fclaw3d_clawpatch46_fort_local_ghost_pack



