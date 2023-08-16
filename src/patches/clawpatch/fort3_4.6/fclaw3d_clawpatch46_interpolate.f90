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

    implicit none
    integer :: mx,my,mz,mbc,meqn,refratio,igrid,idir,iface_coarse
    integer :: num_neighbors
    integer*8 :: transform_ptr
    double precision ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq,r2, m, k
    integer :: ibc, i1
    integer :: jbc, j1
    integer :: ic, jc, mth
    double precision :: gradx, grady, qc, sl, sr, value
    double precision :: fclaw2d_clawpatch_compute_slopes

    !! # This should be refratio*refratio.
    integer :: rr2
    parameter(rr2 = 4)
    integer :: i2(0:rr2-1),j2(0:rr2-1)
    logical :: fclaw2d_clawpatch_is_valid_interp
    logical :: skip_this_grid

    integer :: a(2,2), f(2)
    integer :: ii,jj,dc(2),df(2,0:rr2-1),iff,jff
    double precision :: shiftx(0:rr2-1),shifty(0:rr2-1)

    !! exit with error
    STOP 'NOT IMPLIMENTED'

    mth = 5
    r2 = refratio*refratio
    if (r2 .ne. rr2) then
        write(6,*) 'average_face_ghost (claw2d_utils.f) ', & 
        '  Refratio**2 is not equal to rr2'
        stop
    endif

    call fclaw3d_clawpatch_build_transform(transform_ptr,a,f)

    !! # This needs to be written for refratios .ne. 2.
    m = 0
    do jj = 0,1
        do ii = 0,1
            !! # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

            !! # Direction on fine grid (converted using metric). Divide
            !! # by refratio to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/refratio
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/refratio

            !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5d0)/2.d0
            shifty(m) = (jj-0.5d0)/2.d0
            m = m + 1
        enddo
    enddo
    !! # Create map :

    mq_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            if (idir .eq. 0) then
                !! # this ensures that we get 'hanging' corners        
                do ibc = 1,mbc/2
                    if (iface_coarse .eq. 0) then
                        ic = ibc
                    elseif (iface_coarse .eq. 1) then
                        ic = mx - ibc + 1
                    else
                        write(6,*) 'interpolate : Problem with iface_coarse'
                        write(6,*) 'iface_coarse = ', iface_coarse
                        stop               
                    endif
                    do jc = 1,mx
                        i1 = ic
                        j1 = jc
                        call fclaw3d_clawpatch_transform_face_half(i1,j1,i2,j2,transform_ptr)
                        skip_this_grid = .false.
                        do m = 0,r2-1
                            if (.not. fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc)) then
                                skip_this_grid = .true.
                                exit
                            endif
                        enddo
                        if (.not. skip_this_grid) then
                            qc = qcoarse(ic,jc,k, mq)
                            !! # Compute limited slopes in both x and y. Note we are not
                            !! # really computing slopes, but rather just differences.
                            !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                            sl = (qc - qcoarse(ic-1,jc,k,mq))
                            sr = (qcoarse(ic+1,jc,k,mq) - qc)
                            gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            sl = (qc - qcoarse(ic,jc-1,k,mq))
                            sr = (qcoarse(ic,jc+1,k,mq) - qc)
                            grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            do m = 0,rr2-1
                                iff = i2(0) + df(1,m)
                                jff = j2(0) + df(2,m)
                                value = qc + gradx*shiftx(m) + grady*shifty(m)
                                qfine(iff,jff,k,mq) = value
                            enddo
                        endif
                    enddo !! jc loop
                enddo !! ibc loop            
            else
                do jbc = 1,mbc/2
                    if (iface_coarse .eq. 2) then
                        jc = jbc
                    elseif (iface_coarse .eq. 3) then
                        !! # iface_coarse = 3
                        jc = my - jbc + 1
                    else
                        write(6,*) 'interpolate : Problem with iface_coarse'
                        write(6,*) 'iface_coarse = ', iface_coarse
                        stop
                    endif
                    do ic = 1,mx
                        i1 = ic
                        j1 = jc
                        call fclaw3d_clawpatch_transform_face_half(i1,j1,i2,j2, transform_ptr)
                        !! # ---------------------------------------------
                        !! # Two 'half-size' neighbors will be passed into
                        !! # this routine.  Only half of the coarse grid ghost
                        !! # indices will be valid for the particular grid
                        !! # passed in.  We skip those ghost cells that will
                        !! # have to be filled in by the other half-size
                        !! # grid.
                        !! # ---------------------------------------------
                        skip_this_grid = .false.
                        do m = 0,r2-1
                            if (.not. fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc)) then
                               skip_this_grid = .true.
                               exit
                            endif
                        enddo
                        if (.not. skip_this_grid) then
                            qc = qcoarse(ic,jc,k,mq)

                            sl = (qc - qcoarse(ic-1,jc,k,mq))
                            sr = (qcoarse(ic+1,jc,k,mq) - qc)
                            gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            sl = (qc - qcoarse(ic,jc-1,k,mq))
                            sr = (qcoarse(ic,jc+1,k,mq) - qc)
                            grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                            do m = 0,rr2-1
                                iff = i2(0) + df(1,m)
                                jff = j2(0) + df(2,m)
                                value = qc + gradx*shiftx(m) + grady*shifty(m)
                                qfine(iff,jff,k,mq) = value                                
                            enddo
                        endif                    
                    enddo  !! ic loop
                enddo !! jbc loop
            endif !! end idir branch
        enddo k_loop
    enddo mq_loop

end subroutine fclaw3d_clawpatch46_fort_interpolate_face

subroutine fclaw3d_clawpatch46_fort_interpolate_corner & 
           (mx,my,mz, mbc,meqn,refratio, & 
           qcoarse,qfine,icorner_coarse,transform_ptr)
    implicit none

    integer mx,my,mz,mbc,meqn,icorner_coarse,refratio
    integer*8 transform_ptr
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer ic, jc, mq, ibc,jbc, mth, k
    double precision qc, sl, sr, gradx, grady
    double precision fclaw2d_clawpatch_compute_slopes, value

    !! # This should be refratio*refratio.
    integer i1,j1,m, r2
    integer rr2
    parameter(rr2 = 4)
    integer i2(0:rr2-1),j2(0:rr2-1)

    integer a(2,2), f(2)
    integer ii,jj,iff,jff,dc(2),df(2,0:rr2-1)
    double precision shiftx(0:rr2-1), shifty(0:rr2-1)

    !! exit with error
    STOP 'NOT IMPLIMENTED'

    r2 = refratio*refratio
    if (r2 .ne. rr2) then
        write(6,*) 'average_corner_ghost (claw2d_utils.f) ', & 
                  '  Refratio**2 is not equal to rr2'
        stop
    endif

    call fclaw3d_clawpatch_build_transform(transform_ptr,a,f)

    m = 0
    do jj = 0,1
        do ii = 0,1
            !! # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

            !! # Direction on fine grid (converted using metric). Divide
            !! # by 2 (refratio) to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/2
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/2

            !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5d0)/2.d0
            shifty(m) = (jj-0.5d0)/2.d0
            m = m + 1
        enddo
    enddo


    mth = 5

    do ibc = 1,mbc/2
        do jbc = 1,mbc/2
            if (icorner_coarse .eq. 0) then
                ic = ibc
                jc = jbc
            elseif (icorner_coarse .eq. 1) then
                ic = mx - ibc + 1
                jc = jbc
            elseif (icorner_coarse .eq. 2) then
                ic = ibc
                jc = my - jbc + 1
            elseif (icorner_coarse .eq. 3) then
                ic = mx - ibc + 1
                jc = my - jbc + 1
            else
                write(6,*) "interpolate : Problem with icorner_coarse"
                write(6,*) "icorner_coarse = ", icorner_coarse
                stop                
            endif

            !! # Interpolate coarse grid corners to fine grid corner ghost cells
            i1 = ic
            j1 = jc
            call fclaw3d_clawpatch_transform_corner_half(i1,j1,i2,j2, transform_ptr)

            mq_loop : do mq = 1,meqn
                k_loop : do k = 1,mz
                    qc = qcoarse(ic,jc,k,mq)

                    !! # Compute limited slopes in both x and y. Note we are not
                    !! # really computing slopes, but rather just differences.
                    !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                    sl = (qc - qcoarse(ic-1,jc,k,mq))
                    sr = (qcoarse(ic+1,jc,k,mq) - qc)
                    gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc-1,k,mq))
                    sr = (qcoarse(ic,jc+1,k,mq) - qc)
                    grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    do m = 0,rr2-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m)
                        qfine(iff,jff,k,mq) = value
                    end do
                end do k_loop
            end do mq_loop
        end do !! jbc_loop
    end do !! ibc_loop

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

    integer :: ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth, k
    integer :: ic,jc,ic_add, jc_add, ifine, jfine
    double precision :: qc, shiftx, shifty, sl, sr, gradx, grady
    double precision :: fclaw2d_clawpatch_compute_slopes

    integer :: p4est_refineFactor,refratio

    !! exit with error
    STOP 'NOT IMPLIMENTED'

    p4est_refineFactor = 2
    refratio = 2

    !! # Use limiting done in AMRClaw.
    mth = 5

    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    i1 = 1-ig
    i2 = mx/p4est_refineFactor + (1-ig)
    ic_add = ig*mx/p4est_refineFactor

    j1 = 1-jg
    j2 = my/p4est_refineFactor + (1-jg)
    jc_add = jg*my/p4est_refineFactor

    mq_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            do j = j1,j2
                do i = i1,i2
                    ic = i + ic_add
                    jc = j + jc_add
                    qc = qcoarse(ic,jc,k,mq)

                    !! # Compute limited slopes in both x and y. Note we are not
                    !! # really computing slopes, but rather just differences.
                    !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                    sl = (qc - qcoarse(ic-1,jc,k,mq))
                    sr = (qcoarse(ic+1,jc,k,mq) - qc)
                    gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    sl = (qc - qcoarse(ic,jc-1,k,mq))
                    sr = (qcoarse(ic,jc+1,k,mq) - qc)
                    grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                    !! # Fill in refined values on coarse grid cell (ic,jc)
                    do ii = 1,refratio
                        do jj = 1,refratio
                            shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                            shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                            ifine = (i-1)*refratio + ii
                            jfine = (j-1)*refratio + jj  
                            qfine(ifine,jfine,k,mq) = qc + shiftx*gradx + shifty*grady
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


