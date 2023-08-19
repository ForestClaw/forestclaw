!! # ----------------------------------------------------------
!! # Interpolation routines - (i,j,mq) ordering
!! # ----------------------------------------------------------
!! # interpolate_face_ghost
!! # interpolate_corner_ghost
!! # interpolate_to_fine_patch
!! #
!! # Other routines :
!! # fclaw2d_clawpatch_compute_slopes 
!! # (for limited function reconstruction)
!! # fixcapaq (to preserve conservation)
!! #
!! # Note that fixcapaq is only used when regridding;  ghost
!! # cell interpolation is not conservative in the mapped case.
!! # (Should it be?  We are going to correct the flux mixmatch
!! # anyhow, so maybe the accuracy of the ghost cell values is
!! # more important.)
!! # ----------------------------------------------------------



!! # ----------------------------------------------------------
!! # This routine is used for both mapped and non-mapped
!! # cases.
!! # ----------------------------------------------------------
subroutine fclaw3d_clawpatch46_fort_interpolate_face & 
    (mx,my,mz,mbc,meqn,qcoarse,qfine, & 
     idir,iface_coarse,num_neighbors,refratio,igrid, & 
     transform_ptr)
    use fclaw3d_clawpatch_utils

    implicit none
    integer :: mx,my,mz,mbc,meqn,refratio,igrid,idir,iface_coarse
    integer :: num_neighbors
    integer*8 :: transform_ptr
    double precision ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq,r3, m
    integer :: ic, jc, kc,  mth
    integer :: i_start, j_start, k_start, i_end, j_end, k_end
    double precision :: gradx, grady, gradz, qc, sl, sr, value
    double precision :: fclaw2d_clawpatch_compute_slopes

    !! # This should be refratio*refratio.
    integer :: rr3
    parameter(rr3 = 8)
    integer :: i2(0:rr3-1),j2(0:rr3-1),k2(0:rr3-1)
    logical :: fclaw3d_clawpatch_is_valid_interp
    logical :: skip_this_grid

    integer :: a(3,3), f(3)
    integer :: ii,jj,kk,dc(3),df(3,0:rr3-1),iff,jff,kff
    double precision :: shiftx(0:rr3-1),shifty(0:rr3-1),shiftz(0:rr3-1)

    mth = 5
    r3 = refratio**3
    if (r3 .ne. rr3) then
        write(6,*) 'average_face_ghost (claw2d_utils.f) ', & 
        '  Refratio**2 is not equal to rr2'
        stop
    endif

    call build_transform(transform_ptr,a,f)

    !! # This needs to be written for refratios .ne. 2.
    m = 0
    do kk = 0,1
        do jj = 0,1
            do ii = 0,1
                !! # Direction on coarse grid
                dc(1) = ii
                dc(2) = jj
                dc(3) = kk

                !! # Direction on fine grid (converted using metric). Divide
                !! # by refratio to scale length to unit vector
                df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2) + a(1,3)*dc(3))/refratio
                df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2) + a(2,3)*dc(3))/refratio
                df(3,m) = (a(3,1)*dc(1) + a(3,2)*dc(2) + a(3,3)*dc(3))/refratio

                !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
                shiftx(m) = (ii-0.5d0)/2.d0
                shifty(m) = (jj-0.5d0)/2.d0
                shiftz(m) = (kk-0.5d0)/2.d0
                m = m + 1
            enddo
        enddo
    enddo
    !! # Create map :

    !! loop over coarse iface interior
    call interior_face_bounds( &
        iface_coarse,mx,my,mz,mbc, &
        i_start, j_start, k_start, &
        i_end,   j_end,   k_end)

    mq_loop : do mq = 1,meqn
        k_loop : do kc = k_start,k_end
            j_loop : do jc = j_start,j_end
                i_loop : do ic = i_start,i_end 
                        call fclaw3d_clawpatch_transform_face_half(ic,jc,kc,i2,j2,k2,transform_ptr)
                        skip_this_grid = .false.
                        do m = 0,r3-1
                            if (.not. is_valid_interp(i2(m),j2(m),k2(m),mx,my,mz,mbc)) then
                                skip_this_grid = .true.
                                exit
                            endif
                        enddo
                        if (.not. skip_this_grid) then
                            qc = qcoarse(ic,jc,kc, mq)
                            !! # Compute limited slopes in both x and y. Note we are not
                            !! # really computing slopes, but rather just differences.
                            !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                            sl = (qc - qcoarse(ic-1,jc,kc,mq))
                            sr = (qcoarse(ic+1,jc,kc,mq) - qc)
                            gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            sl = (qc - qcoarse(ic,jc-1,kc,mq))
                            sr = (qcoarse(ic,jc+1,kc,mq) - qc)
                            grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            sl = (qc - qcoarse(ic,jc,kc-1,mq))
                            sr = (qcoarse(ic,jc,kc+1,mq) - qc)
                            gradz = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            do m = 0,rr3-1
                                iff = i2(0) + df(1,m)
                                jff = j2(0) + df(2,m)
                                kff = k2(0) + df(3,m)
                                value = qc + gradx*shiftx(m) + grady*shifty(m) + gradz*shiftz(m)
                                qfine(iff,jff,kff,mq) = value
                            enddo
                        endif
                enddo i_loop
            enddo j_loop
        enddo k_loop
    enddo mq_loop

end subroutine fclaw3d_clawpatch46_fort_interpolate_face

subroutine fclaw3d_clawpatch46_fort_interpolate_corner & 
           (mx,my,mz, mbc,meqn,refratio, & 
           qcoarse,qfine,icorner_coarse,transform_ptr)
    use fclaw3d_clawpatch_utils
    implicit none

    integer mx,my,mz,mbc,meqn,icorner_coarse,refratio
    integer*8 transform_ptr
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer mq, ic,jc,kc, mth
    double precision qc, sl, sr, gradx, grady, gradz
    double precision fclaw2d_clawpatch_compute_slopes, value
    integer i_start, j_start, k_start

    !! # This should be refratio*refratio.
    integer m, r3
    integer rr3
    parameter(rr3 = 8)
    integer i2(0:rr3-1),j2(0:rr3-1),k2(0:rr3-1)

    integer a(3,3), f(3)
    integer ii,jj,kk,iff,jff,kff,dc(3),df(3,0:rr3-1)
    double precision shiftx(0:rr3-1), shifty(0:rr3-1), shiftz(0:rr3-1)

    r3 = refratio**3
    if (r3 .ne. rr3) then
        write(6,*) 'average_corner_ghost (claw2d_utils.f) ', & 
                  '  Refratio**3 is not equal to rr3'
        stop
    endif

    call build_transform(transform_ptr,a,f)

    m = 0
    do kk = 0,1
        do jj = 0,1
            do ii = 0,1
                !! # Direction on coarse grid
                dc(1) = ii
                dc(2) = jj
                dc(3) = kk

                !! # Direction on fine grid (converted using metric). Divide
                !! # by 2 (refratio) to scale length to unit vector
                df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2) + a(1,3)*dc(3))/2
                df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2) + a(2,3)*dc(3))/2
                df(3,m) = (a(3,1)*dc(1) + a(3,2)*dc(2) + a(3,3)*dc(3))/2

                !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
                shiftx(m) = (ii-0.5d0)/2.d0
                shifty(m) = (jj-0.5d0)/2.d0
                shiftz(m) = (kk-0.5d0)/2.d0
                m = m + 1
            enddo
        enddo
    enddo


    mth = 5

    !! get lower-bottom-left corner of coarse ghost cells that 
    !! we are interpolating from
    call fclaw3d_clawpatch46_fort_get_corner_start_coarse_to_fine( &
            icorner_coarse,refratio, &
            mx,my,mz,mbc, &
            i_start,j_start,k_start)

    mq_loop : do mq = 1,meqn
        kc_loop : do kc = k_start,k_start+mbc/2-1
            jc_loop : do jc = j_start,j_start+mbc/2-1
                ic_loop : do ic = i_start,i_start+mbc/2-1
                    !! # Interpolate coarse grid corners to fine grid corner ghost cells
                    call fclaw3d_clawpatch_transform_corner_half(ic,jc,kc,i2,j2,k2, transform_ptr)

                    qc = qcoarse(ic,jc,kc,mq)

                    !! # Compute limited slopes in both x and y. Note we are not
                    !! # really computing slopes, but rather just differences.
                    !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                    sl = (qc - qcoarse(ic-1,jc,kc,mq))
                    sr = (qcoarse(ic+1,jc,kc,mq) - qc)
                    gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc-1,kc,mq))
                    sr = (qcoarse(ic,jc+1,kc,mq) - qc)
                    grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc,kc-1,mq))
                    sr = (qcoarse(ic,jc,kc+1,mq) - qc)
                    gradz = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    do m = 0,rr3-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        kff = k2(0) + df(3,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m) + gradz*shiftz(m)
                        qfine(iff,jff,kff,mq) = value
                    end do
                end do ic_loop
            end do jc_loop
        end do kc_loop
    end do mq_loop

end subroutine fclaw3d_clawpatch46_fort_interpolate_corner


!! # Conservative intepolation to fine grid patch
subroutine fclaw3d_clawpatch46_fort_interpolate2fine & 
          (mx,my,mz,mbc,meqn,qcoarse, qfine, volcoarse, &
           volfine, igrid, manifold)
    implicit none

    integer :: mx,my,mz,mbc,meqn
    integer :: igrid, manifold

    double precision :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    double precision :: volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision ::   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    integer :: ii, jj, kk, i,j,k, i1, i2, j1, j2, k1, k2, ig, jg, kg, mq, mth
    integer :: ic,jc,kc,ic_add, jc_add, kc_add, ifine, jfine, kfine
    integer :: i_start_fine, j_start_fine, k_start_fine
    double precision :: qc, shiftx, shifty, shiftz, sl, sr, gradx, grady, gradz
    double precision :: fclaw2d_clawpatch_compute_slopes

    logical :: lower_x, lower_y, lower_z

    integer :: p8est_refineFactor,refratio

    p8est_refineFactor = 2
    refratio = 2

    !! # Use limiting done in AMRClaw.
    mth = 5

    ig = merge(1,0,btest(igrid,0))
    jg = merge(1,0,btest(igrid,1))
    kg = merge(1,0,btest(igrid,2))

    i1 = 1-ig
    i2 = mx/p8est_refineFactor + (1-ig)
    ic_add = ig*mx/p8est_refineFactor

    j1 = 1-jg
    j2 = my/p8est_refineFactor + (1-jg)
    jc_add = jg*my/p8est_refineFactor

    k1 = 1-kg
    k2 = mz/p8est_refineFactor + (1-kg)
    kc_add = kg*mz/p8est_refineFactor

    !!iterate of part of coarse patch that we are intepolating from
    mq_loop : do mq = 1,meqn
        k_loop : do k = k1,k2
            do j = j1,j2
                do i = i1,i2
                    ic = i + ic_add
                    jc = j + jc_add
                    kc = k + kc_add
                    qc = qcoarse(ic,jc,kc,mq)

                    !! # Compute limited slopes in both x and y. Note we are not
                    !! # really computing slopes, but rather just differences.
                    !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                    sl = (qc - qcoarse(ic-1,jc,kc,mq))
                    sr = (qcoarse(ic+1,jc,kc,mq) - qc)
                    gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc-1,kc,mq))
                    sr = (qcoarse(ic,jc+1,kc,mq) - qc)
                    grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc,kc-1,mq))
                    sr = (qcoarse(ic,jc,kc+1,mq) - qc)
                    gradz = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    !! # Fill in refined values on coarse grid cell (ic,jc,kc)
                    do kk = 1,refratio
                        do jj = 1,refratio
                            do ii = 1,refratio
                                shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                                shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                                shiftz = (kk - refratio/2.d0 - 0.5d0)/refratio
                                ifine = (i-1)*refratio + ii
                                jfine = (j-1)*refratio + jj  
                                kfine = (k-1)*refratio + kk  
                                qfine(ifine,jfine,kfine,mq) = qc + shiftx*gradx + shifty*grady + shiftz*gradz
                            end do
                        end do
                    end do
                end do  !! i_loop
            end do !! j_loop
        end do k_loop
    end do mq_loop

    if (manifold .ne. 0) then
        !!write(6,*) 'interpolate:fixcapaq2 : Manifold not yet implemented in 3D'
        !!stop
        call fclaw3d_clawpatch46_fort_fixcapaq2(mx,my,mz,mbc,meqn, & 
                        qcoarse,qfine, volcoarse,volfine,igrid)
    endif


end subroutine  fclaw3d_clawpatch46_fort_interpolate2fine


    !! # ------------------------------------------------------
    !! # So far, this is only used by the interpolation from
    !! # coarse to fine when regridding.  But maybe it should
    !! # be used by the ghost cell routines as well?
    !! # ------------------------------------------------------
subroutine fclaw3d_clawpatch46_fort_fixcapaq2(mx,my,mz,mbc,meqn, & 
           qcoarse,qfine, volcoarse,volfine,igrid)
    implicit none

    integer :: mx,my,mz,mbc,meqn, refratio, igrid
    integer :: p4est_refineFactor

    double precision ::  qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::    qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision ::    volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    integer :: i,j,k,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
    double precision :: kf, kc, r2, sum, cons_diff, qf, qc, volf, dz, volc

    p4est_refineFactor = 2
    refratio = 2

    !! This is a bogus value, since we shouldn't end up here (yet). 
    dz = 1

    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/p4est_refineFactor
    jc_add = jg*my/p4est_refineFactor

    !! # ------------------------------------------------------
    !! # This routine ensures that the interpolated solution
    !! # has the same mass as the coarse grid solution
    !! # -------------------------------------------------------

    r2 = refratio*refratio
    mq_loop : do m = 1,meqn
        k_loop : do k = 1,mz
            do i = 1,mx/p4est_refineFactor
                do j = 1,my/p4est_refineFactor
                    sum = 0.d0
                    do ii = 1,refratio
                        do jj = 1,refratio
                           ifine = (i-1)*refratio + ii
                           jfine = (j-1)*refratio + jj
                           kf = volfine(ifine,jfine,k)
                           volf = kf*dz
                           qf = qfine(ifine,jfine,k,m)
                           sum = sum + volf*qf
                        enddo
                    enddo

                    kc = volcoarse(i+ic_add,j+jc_add,k)
                    volc = kc*dz
                    qc = qcoarse(i+ic_add, j+jc_add,k,m)
                    cons_diff = (qc*volc - sum)/r2

                    do ii = 1,refratio
                        do jj = 1,refratio
                           ifine  = (i-1)*refratio + ii
                           jfine  = (j-1)*refratio + jj
                           kf = volfine(ifine,jfine,k)
                           volf = kf*dz
                           qfine(ifine,jfine,k,m) = qfine(ifine,jfine,k,m) + cons_diff/volf
                       end do
                    end do
                end do !! j loop
            enddo !! i loop
        enddo k_loop
    enddo mq_loop

end


