!!> \file
!!> \defgroup Averaging Average fine grids to a coarse grid
!!> Average cells from coarse grid to fine grid.
!!>
!!> Routines described here average are used to fill coarse grid ghost
!!> cells, and average sibling grids onto a parent grid.  Indices
!!> for cells at block boundaries are transformed using encodings
!!> stored in `transform_cptr`.
!!>
!!> \param [in] mx,my       Number of cells in x,y direction
!!> \param [in] mbc      Number of ghost cells
!!> \param [in] meqn     Number of equations
!!> \param [in] qcoarse,qfine  Solution on coarse,fine grid
!!> \param [in] areacoarse,area  Area of mesh cells on coarse,fine grids.
!!> \param [in] idir     Face orientation - 0 for x-faces; 1 for y-faces [0-1]
!!> \param [in] iface    Face number of fine grid [0-3].
!!> \param [in] iface_coarse Face number of coarse grid [0-3].
!!> \param [in] num_neighbors Number of fine grid neighbors [2].
!!> \param [in] refratio  Refinement ratio between coarse and fine grids [2].
!!> \param [in] manifold  Flag indicating whether we are on mapped grid [0-1].
!!> \param [in] transform_cptr  Encoding for indices at block boundaries (C only).

!!> \ingroup Averaging
!!> Average fine ghost cell values.
!!>
!!> Average fine grid interior values to neighboring ghost cell values of
!!> the coarse grid.

subroutine fclaw3d_clawpatch_face_ghost_bounds( &
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
            i_start = mx+1
            i_end = mx+mbc
        else
            i_start = 1-mbc
            i_end = 0
        endif
        j_start = 1
        j_end = my
        k_start = 1
        k_end = mz
    else if(axis .eq. 1) then
        i_start = 0
        i_end = mx
        if (upper) then
            j_start = my+1
            j_end = my+mbc
        else
            j_start = 1-mbc
            j_end = 0
        endif
        k_start = 1
        k_end = mz
    else if(axis .eq. 2) then
        i_start = 1
        i_end = mx
        j_start = 1
        j_end = my
        if (upper) then
            k_start = mz+1
            k_end = mz+mbc
        else
            k_start = 1-mbc
            k_end = 0
        endif
    endif
end subroutine fclaw3d_clawpatch_face_ghost_bounds

    
SUBROUTINE fclaw3d_clawpatch46_fort_average_face(mx,my,mz,mbc,meqn, & 
           qcoarse,qfine,volcoarse, volfine, & 
           idir,iface_coarse,num_neighbors,refratio,igrid, & 
           manifold, transform_cptr)
    USE fclaw3d_clawpatch_utils
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn,refratio,igrid,idir,iface_coarse
    INTEGER :: manifold
    INTEGER*8 :: transform_cptr
    INTEGER :: num_neighbors
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold.
    DOUBLE PRECISION :: volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    DOUBLE PRECISION ::   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    INTEGER :: mq,r3, m
    INTEGER :: ic, jc, kc
    INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end

    !! # For the extruded mesh, this is still only rr2
    INTEGER :: rr3
    PARAMETER(rr3 = 8)
    INTEGER, DIMENSION(0:rr3-1) :: i2, j2, k2
    !! DOUBLE PRECISION :: kc

    LOGICAL :: skip_this_grid
    DOUBLE PRECISION :: vf_sum
    DOUBLE PRECISION :: sum, qf, kf
    LOGICAL :: is_manifold


    is_manifold = manifold .eq. 1

    r3 = refratio**3    
    if (r3 .ne. rr3) then
        write(6,*) 'average_face_ghost (fclaw3d_clawpatch46_average.f90) ', & 
             '  Refratio**3 is not equal to rr3'
        stop
    endif

    !! # Average fine grid onto coarse grid
    call fclaw3d_clawpatch_face_ghost_bounds( &
        iface_coarse,mx,my,mz,mbc, &
        i_start, j_start, k_start, &
        i_end,   j_end,   k_end)

    meqn_loop : do mq = 1,meqn
        k_loop : do kc = k_start,k_end
            j_loop : do jc = j_start,j_end
                i_loop : do ic = i_start,i_end
                        call fclaw3d_clawpatch_transform_face_half(ic,jc,kc,i2,j2,k2,transform_cptr)
                        !! # ---------------------------------------------
                        !! # Two 'half-size' neighbors will be passed into
                        !! # this routine.  Only half of the coarse grid ghost
                        !! # indices will be valid for the particular grid
                        !! # passed in.  We skip those ghost cells that will
                        !! # have to be filled in by the other half-size
                        !! # grid.
                        !! # ---------------------------------------------
                        skip_this_grid = .false.
                        do m = 0,r3-1
                            if (.not. is_valid_average(i2(m),j2(m),k2(m),mx,my,mz)) then
                                skip_this_grid = .true.
                                exit 
                            endif
                        end do
                        if (.not. skip_this_grid) then
                            if (is_manifold) then
                                sum = 0
                                vf_sum = 0
                                do m = 0,r3-1
                                    qf = qfine(i2(m),j2(m),k2(m),mq)
                                    kf = volfine(i2(m),j2(m),k2(m))
                                    sum = sum + qf*kf
                                    vf_sum = vf_sum + kf
                                end do
                                !! # Use vols of the fine grid mesh cells instead.
                                qcoarse(ic,jc,kc,mq) = sum/vf_sum
                            else
                                sum = 0
                                do m = 0,r3-1
                                    sum = sum + qfine(i2(m),j2(m),k2(m),mq)
                                end do
                                qcoarse(ic,jc,kc,mq) = sum/dble(r3)                            
                            endif
                        endif
                enddo i_loop
            end do j_loop
        end do k_loop
    end do meqn_loop
end subroutine  fclaw3d_clawpatch46_fort_average_face

!!> \ingroup Averaging
!!> Average across edge.

subroutine fclaw3d_clawpatch46_fort_average_edge(mx,my,mz,mbc,meqn, &
    refratio,qcoarse,qfine,volcoarse,volfine, & 
    manifold,iedge_coarse,transform_cptr)
    USE fclaw3d_clawpatch_utils
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn,refratio,iedge_coarse, manifold
    INTEGER*8 :: transform_cptr
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold.
    DOUBLE PRECISION :: volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    DOUBLE PRECISION ::   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    INTEGER :: ibc,jbc,kbc,mq,r3
    LOGICAL :: is_manifold
    DOUBLE PRECISION :: qf,kf, sum

    !! # This should be refratio*refratio.
    INTEGER :: rr3
    PARAMETER(rr3 = 8)
    INTEGER :: i2(0:rr3-1),j2(0:rr3-1),k2(0:rr3-1)

    INTEGER :: m
    INTEGER :: i_start, j_start, k_start, i_end, j_end, k_end
    DOUBLE PRECISION :: vf_sum

    LOGICAL :: skip_this_grid


    r3 = refratio**3
    if (r3 .ne. rr3) then
        write(6,*) 'average_corner_ghost (claw2d_utils.f) ', & 
        '  Refratio**3 is not equal to rr3'
        stop
    endif

    is_manifold = manifold .eq. 1

    !! get lower-bottom-left corner of coarse ghost cells
    call fclaw3d_clawpatch46_fort_get_edge_bounds(iedge_coarse, &
                mx,my,mz,mbc, &
                i_start,j_start,k_start, &
                i_end,  j_end,  k_end)

    meqn_loop : do mq = 1,meqn
        kbc_loop : do kbc = k_start,k_end
            jbc_loop : do jbc = j_start,j_end
                ibc_loop : do ibc = i_start,i_end
                    !! # Average fine grid corners onto coarse grid ghost corners
                    call fclaw3d_clawpatch_transform_edge_half(ibc,jbc,kbc,i2,j2,k2, transform_cptr)
                    skip_this_grid = .false.
                    do m = 0,r3-1
                        if (.not. is_valid_average(i2(m),j2(m),k2(m),mx,my,mz)) then
                            skip_this_grid = .true.
                            exit 
                        endif
                    end do
                    if (.not. skip_this_grid) then
                        if (is_manifold) then
                            sum = 0
                            vf_sum = 0
                            do m = 0,r3-1
                                qf = qfine(i2(m),j2(m),k2(m),mq)
                                kf = volfine(i2(m),j2(m),k2(m))
                                sum = sum + kf*qf
                                vf_sum = vf_sum + kf
                            enddo
                            qcoarse(ibc,jbc,kbc,mq) = sum/vf_sum
                        else
                            sum = 0
                            do m = 0,r3-1
                                qf = qfine(i2(m),j2(m),k2(m),mq)
                                sum = sum + qf
                            end do
                            qcoarse(ibc,jbc,kbc,mq) = sum/dble(r3)
                        endif
                    endif
                enddo ibc_loop
            enddo jbc_loop
        enddo kbc_loop
    end do meqn_loop

end subroutine fclaw3d_clawpatch46_fort_average_edge



!!> \ingroup Averaging
!!> Average across corners.

subroutine fclaw3d_clawpatch46_fort_average_corner(mx,my,mz,mbc,meqn, &
    refratio,qcoarse,qfine,volcoarse,volfine, & 
    manifold,icorner_coarse,transform_cptr)
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn,refratio,icorner_coarse, manifold
    INTEGER*8 :: transform_cptr
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold.
    DOUBLE PRECISION :: volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    DOUBLE PRECISION ::   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    INTEGER :: ibc,jbc,kbc,mq,r3
    LOGICAL :: is_manifold
    DOUBLE PRECISION :: qf,kf, sum

    !! # This should be refratio*refratio.
    INTEGER :: rr3
    PARAMETER(rr3 = 8)
    INTEGER :: i2(0:rr3-1),j2(0:rr3-1),k2(0:rr3-1)

    INTEGER :: m
    INTEGER :: i_start, j_start, k_start
    DOUBLE PRECISION :: vf_sum


    r3 = refratio**3
    if (r3 .ne. rr3) then
        write(6,*) 'average_corner_ghost (claw2d_utils.f) ', & 
        '  Refratio**3 is not equal to rr3'
        stop
    endif

    is_manifold = manifold .eq. 1

    !! get lower-bottom-left corner of coarse ghost cells
    call fclaw3d_clawpatch46_fort_get_corner_start(icorner_coarse, &
                mx,my,mz,mbc, &
                i_start,j_start,k_start)

    meqn_loop : do mq = 1,meqn
        kbc_loop : do kbc = k_start,k_start+mbc-1
            jbc_loop : do jbc = j_start,j_start+mbc-1
                ibc_loop : do ibc = i_start,i_start+mbc-1
                    !! # Average fine grid corners onto coarse grid ghost corners
                    call fclaw3d_clawpatch_transform_corner_half(ibc,jbc,kbc,i2,j2,k2, transform_cptr)
                    if (is_manifold) then
                        sum = 0
                        vf_sum = 0
                        do m = 0,r3-1
                            qf = qfine(i2(m),j2(m),k2(m),mq)
                            kf = volfine(i2(m),j2(m),k2(m))
                            sum = sum + kf*qf
                            vf_sum = vf_sum + kf
                        enddo
                        qcoarse(ibc,jbc,kbc,mq) = sum/vf_sum
                    else
                        sum = 0
                        do m = 0,r3-1
                            qf = qfine(i2(m),j2(m),k2(m),mq)
                            sum = sum + qf
                        end do
                        qcoarse(ibc,jbc,kbc,mq) = sum/dble(r3)
                    endif
                enddo ibc_loop
            enddo jbc_loop
        enddo kbc_loop
    end do meqn_loop

end subroutine fclaw3d_clawpatch46_fort_average_corner


!!> \ingroup  Averaging
!!> Average fine grid siblings to parent coarse grid.
subroutine fclaw3d_clawpatch46_fort_average2coarse(mx,my,mz,mbc,meqn, & 
           qcoarse,qfine,volcoarse,volfine, igrid,manifold)
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn
    INTEGER :: manifold
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold, and so shouldn't
    !! # be referenced. 
    DOUBLE PRECISION :: volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    DOUBLE PRECISION ::   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    !! # This should be refratio**3
    INTEGER :: rr3
    PARAMETER(rr3 = 8)
    INTEGER :: ifine(0:rr3-1),jfine(0:rr3-1),kfine(0:rr3-1)


    INTEGER :: i,j,k,mq, ii,jj,kk, ic, jc, kc, r3, m
    INTEGER :: ig, jg, kg, ic_add, jc_add, kc_add
    INTEGER :: refratio, igrid
    LOGICAL :: is_manifold
    DOUBLE PRECISION :: sum, kf, qf, vf_sum

    refratio = 2

    is_manifold = manifold .eq. 1

    !! # 'iface' is relative to the coarse grid

    r3 = refratio**3
    if (r3 .ne. rr3) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ', & 
         '  Refratio**2 is not equal to rr2'
        stop
    endif


    !! Get (ig,jg,kg) for grid from linear (igrid) coordinates
    !! igrid = ig + refratio*jg + refratio*refratio*kg
    ig = mod(igrid,refratio)
    jg = mod((igrid-ig)/refratio,refratio)
    kg = (igrid-ig-refratio*jg)/(refratio**2)

    !! Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/refratio
    jc_add = jg*my/refratio
    kc_add = kg*mz/refratio

    r3 = refratio**3
    !! this assumes mz,my,mz are divisible by refratio
    !! should probably add a check to clawpatch_options
    meqn_loop : do mq = 1,meqn
        k_loop : do k = 1,mz/refratio
            j_loop : do j = 1,my/refratio
                i_loop : do i = 1,mx/refratio
                    ic = i+ic_add
                    jc = j+jc_add
                    kc = k+kc_add
                    m = 0
                    do kk = 1,refratio
                        do jj = 1,refratio
                            do ii = 1,refratio
                               ifine(m) = (i-1)*refratio + ii
                               jfine(m) = (j-1)*refratio + jj
                               kfine(m) = (k-1)*refratio + kk
                               m = m + 1
                            enddo
                        enddo
                    enddo
                    print *,''
                    if (is_manifold) then
                        sum = 0
                        vf_sum = 0
                        do m = 0,r3-1
                            qf = qfine(ifine(m),jfine(m),kfine(m),mq)
                            kf = volfine(ifine(m),jfine(m),kfine(m))
                            sum = sum + kf*qf
                            vf_sum = vf_sum + kf
                        enddo
                        !!vc = volcoarse(ic,jc,kc)
                        qcoarse(ic,jc,kc,mq) = sum/vf_sum
                    else
                        sum = 0
                        do m = 0,r3-1
                           qf = qfine(ifine(m),jfine(m),kfine(m),mq)
                           sum = sum + qf
                        enddo
                        qcoarse(ic,jc,kc,mq) = sum/r3
                    endif 
                end do i_loop
            enddo j_loop
        enddo k_loop
    enddo meqn_loop
end subroutine  fclaw3d_clawpatch46_fort_average2coarse


