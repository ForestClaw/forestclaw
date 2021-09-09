!! # ----------------------------------------------------------------------
!! # This handles the boundary conditions at the block
!! # corners for the pillow sphere.
!! # ----------------------------------------------------------------------

subroutine fclaw2d_pillow46_copy_block_corner3 (mx,my,mz,mbc,meqn, & 
    qthis, qneighbor, icorner, iblock)
    implicit none

    integer :: mx, my, mz,mbc, meqn, icorner, iblock
    double precision :: qthis(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq, ibc, jbc,k

    do mq = 1,meqn
        k_loop : do k = 1,mz
            if (icorner .eq. 0) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        qthis(1-ibc,1-jbc,k,mq) = qneighbor(ibc,jbc,k,mq)
                        qneighbor(1-ibc,1-jbc,k,mq) = qthis(ibc,jbc,k,mq)
                    end do
                end do
            else if (icorner .eq. 1) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        qthis(mx+ibc,1-jbc,k,mq) = qneighbor(mx+1-ibc,jbc,k,mq)
                        qneighbor(mx+ibc,1-jbc,k,mq) = qthis(mx+1-ibc,jbc,k,mq)
                    end do
                end do
            elseif (icorner .eq. 2) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        qthis(1-ibc,my+jbc,k,mq) = qneighbor(ibc,my+1-jbc,k,mq)
                        qneighbor(1-ibc,my+jbc,k,mq) = qthis(ibc,my+1-jbc,k,mq)
                    end do
                end do
            elseif (icorner .eq. 3) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        qthis(mx+ibc,my+jbc,k,mq) = qneighbor(mx+1-ibc,my+1-jbc,k,mq)
                        qneighbor(mx+ibc,my+jbc,k,mq) = qthis(mx+1-ibc,my+1-jbc,k,mq)
                    end do
                end do
            end if
        end do k_loop
    end do
end subroutine fclaw2d_pillow46_copy_block_corner3

subroutine fclaw2d_pillow46_average_block_corner3(mx,my,mz,dz,mbc,meqn, & 
           refratio, qcoarse, qfine, areacoarse, areafine, & 
           icorner,iblock)
    implicit none

    integer :: mx, my, mz,mbc, meqn, icorner, iblock, refratio
    double precision :: dz
    double precision :: areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    double precision ::   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    double precision :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq, ibc, jbc, k,ii, jj, ifine, jfine
    double precision :: sum,qf,kf,kc,volf


    do mq = 1,meqn
       k_loop : do k = 1-mbc,mz+mbc
            if (icorner .eq. 0) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        sum = 0
                        do ii = 1,refratio
                            do jj = 1,refratio
                                ifine = (ibc-1)*refratio + ii
                                jfine = (jbc-1)*refratio + jj
                                qf = qfine(ifine,jfine,k,mq)
                                kf = areafine(ifine,jfine)
                                volf = kf*dz
                                sum = sum + volf*qf
                            end do
                        end do
                        kc = areacoarse(1-ibc,1-ibc)*dz
                        qcoarse(1-ibc,1-jbc,k,mq) = sum/kc
                    end do
                end do
            elseif (icorner .eq. 1) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        sum = 0
                        do ii = 1,refratio
                            do jj = 1,refratio
                                ifine = (ibc-1)*refratio + ii
                                jfine = (jbc-1)*refratio + jj
                                qf = qfine(mx+1-ifine,jfine,k,mq)
                                kf = areafine(mx+1-ifine,jfine)
                                volf = kf*dz
                                sum = sum + volf*qf
                            end do
                        end do
                        kc = areacoarse(mx+ibc,1-jbc)*dz
                        qcoarse(mx+ibc,1-jbc,k,mq) = sum/kc
                    end do
                end do
            elseif (icorner .eq. 2) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        sum = 0
                        do ii = 1,refratio
                            do jj = 1,refratio
                                ifine = (ibc-1)*refratio + ii
                                jfine = (jbc-1)*refratio + jj
                                qf = qfine(ifine,my+1-jfine,k,mq)
                                kf = areafine(ifine,my+1-jfine)
                                volf = kf*dz
                                sum = sum + volf*qf
                            end do
                        end do
                        kc = areacoarse(1-ibc,my+jbc)*dz
                        qcoarse(1-ibc,my+jbc,k,mq) = sum/kc
                    end do
                end do
            elseif (icorner .eq. 3) then
                do ibc = 1,mbc
                    do jbc = 1,mbc
                        sum = 0
                        do ii = 1,refratio
                            do jj = 1,refratio
                                ifine = (ibc-1)*refratio + ii
                                jfine = (jbc-1)*refratio + jj
                                qf = qfine(mx+1-ifine,my+1-jfine,k,mq)
                                kf = areafine(mx+1-ifine,my+1-jfine)
                                volf = kf*dz
                                sum = sum + volf*qf
                            end do
                        end do
                        kc = areacoarse(mx+ibc,my+jbc)*dz
                        qcoarse(mx+ibc,my+jbc,k,mq) = sum/kc
                    end do
                end do
            end if
        end do k_loop
    end do 

end subroutine fclaw2d_pillow46_average_block_corner3


subroutine fclaw2d_pillow46_interpolate_block_corner3(mx,my,mz,mbc,meqn, & 
           refratio, qcoarse, qfine, icorner_coarse, iblock)
    implicit none
    integer :: mx, my, mz, mbc, meqn, icorner_coarse, iblock, refratio
    double precision :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq, ibc, jbc
    integer :: ic, jc, k, mth
    double precision :: gradx, grady, shiftx, shifty
    double precision :: sl, sr, qc, value
    double precision :: fclaw3dx_clawpatch_compute_slopes

    mth = 5

    if (icorner_coarse .eq. 0) then
        ic = 1
        jc = 1
    elseif (icorner_coarse .eq. 1) then
        ic = mx
        jc = 1
    elseif (icorner_coarse .eq. 2) then
       ic = 1
       jc = my
    elseif (icorner_coarse .eq. 3) then
        ic = mx
        jc = my
    else
        write(6,*) 'pillow : icorner_coarse has unexpected value'
        write(6,*) 'icorner_coarse : ', icorner_coarse
        stop
    endif

    !! # This may not even matter
    do mq = 1,meqn
        k_loop : do k = 1,mz
            qc = qcoarse(ic,jc,k,mq)
            sl = (qc - qcoarse(ic-1,jc,k,mq))
            sr = (qcoarse(ic+1,jc,k,mq) - qc)
            gradx = fclaw3dx_clawpatch_compute_slopes(sl,sr,mth)

            sl = (qc - qcoarse(ic,jc-1,k,mq))
            sr = (qcoarse(ic,jc+1,k,mq) - qc)
            grady = fclaw3dx_clawpatch_compute_slopes(sl,sr,mth)

            !! # Loop over fine grid ghost cells
            do ibc = 1,mbc
                do jbc = 1,mbc
                    !! # Fill in interpolated values on fine grid cell
                    shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
                    shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

                    value = qc + shiftx*gradx + shifty*grady
                    if (icorner_coarse .eq. 0) then
                        qfine(1-ibc,1-jbc,k,mq) = value
                    elseif (icorner_coarse .eq. 1) then
                        qfine(mx+ibc,1-jbc,k,mq) = value
                    elseif (icorner_coarse .eq. 2) then
                        qfine(1-ibc,my+jbc,k,mq) = value
                    else if (icorner_coarse .eq. 3) then
                        qfine(mx+ibc,my+jbc,k,mq) = value
                    endif
                end do
            end do
        end do k_loop
    end do
end subroutine fclaw2d_pillow46_interpolate_block_corner3
