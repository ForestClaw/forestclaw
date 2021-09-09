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
      
SUBROUTINE fclaw3dx_clawpatch46_fort_average_face(mx,my,mz,mbc,meqn, & 
           qcoarse,qfine,areacoarse, areafine, & 
           idir,iface_coarse,num_neighbors,refratio,igrid, & 
           manifold, transform_cptr)
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn,refratio,igrid,idir,iface_coarse
    INTEGER :: manifold
    INTEGER*8 :: transform_cptr
    INTEGER :: num_neighbors
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold.
    DOUBLE PRECISION :: areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    DOUBLE PRECISION ::   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    INTEGER :: mq,r2, m
    INTEGER :: ic, ibc, jc, jbc, k

    !! # For the extruded mesh, this is still only rr2
    INTEGER :: rr2
    PARAMETER(rr2 = 4)
    INTEGER, DIMENSION(0:rr2-1) :: i2, j2
    !! DOUBLE PRECISION :: kc

    LOGICAL :: fclaw2d_clawpatch_is_valid_average, skip_this_grid
    DOUBLE PRECISION :: vf_sum
    DOUBLE PRECISION :: sum, qf, kf
    LOGICAL :: is_manifold


    is_manifold = manifold .eq. 1

    !! # 'iface' is relative to the coarse grid

    r2 = refratio**2    
    if (r2 .ne. rr2) then
        write(6,*) 'average_face_ghost (fclaw3dx_clawpatch46_average.f90) ', & 
             '  Refratio**2 is not equal to rr2'
        stop
    endif

    !! # Average fine grid onto coarse grid
    meqn_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            if (idir .eq. 0) then
                do jc = 1,my
                    do ibc = 1,mbc
                        !! # ibc = 1 corresponds to first layer of ghost cells, and
                        !! # ibc = 2 corresponds to the second layer

                        if (iface_coarse .eq. 0) then
                            ic = 1-ibc
                        elseif (iface_coarse .eq. 1) then
                            ic = mx+ibc
                        endif
                        call fclaw3dx_clawpatch_transform_face_half(ic,jc,i2,j2,transform_cptr)
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
                            if (.not. fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my)) then
                                skip_this_grid = .true.
                                exit 
                            endif
                        end do
                        if (.not. skip_this_grid) then
                            if (is_manifold) then
                                sum = 0
                                vf_sum = 0
                                do m = 0,r2-1
                                    qf = qfine(i2(m),j2(m),k,mq)
                                    kf = areafine(i2(m),j2(m))
                                    sum = sum + qf*kf
                                    vf_sum = vf_sum + kf
                                end do
                                !! # Use vols of the fine grid mesh cells instead.
                                qcoarse(ic,jc,k,mq) = sum/vf_sum
                            else
                                sum = 0
                                do m = 0,r2-1
                                    sum = sum + qfine(i2(m),j2(m),k,mq)
                                end do
                                qcoarse(ic,jc,k,mq) = sum/dble(r2)                            
                            endif
                        endif
                    enddo
                end do
            else
                !! # idir = 1 (faces 2,3)
                do jbc = 1,mbc
                    do ic = 1,mx
                        if (iface_coarse .eq. 2) then
                            jc = 1-jbc
                        elseif (iface_coarse .eq. 3) then
                            jc = my+jbc
                        endif

                        call fclaw3dx_clawpatch_transform_face_half(ic,jc,i2,j2, transform_cptr)
                        skip_this_grid = .false.
                        do m = 0,r2-1
                            if (.not. fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my)) then
                                skip_this_grid = .true.
                            endif
                        end do
                        if (.not. skip_this_grid) then
                            if (is_manifold) then
                                sum = 0
                                vf_sum = 0
                                do m = 0,r2-1
                                    qf = qfine(i2(m),j2(m),k, mq)
                                    kf = areafine(i2(m),j2(m))
                                    sum = sum + qf*kf
                                    vf_sum = vf_sum + kf
                                end do
                                qcoarse(ic,jc,k,mq) = sum/vf_sum
                            else
                                sum = 0
                                do m = 0,r2-1
                                    sum = sum + qfine(i2(m),j2(m),k,mq)
                                end do
                                qcoarse(ic,jc,k,mq) = sum/dble(r2)
                            endif            
                        endif               
                    enddo  !! i loop
                enddo  !! jbc loop
            endif  !! idir conditonal
        end do k_loop
    end do meqn_loop
end subroutine  fclaw3dx_clawpatch46_fort_average_face


!!> \ingroup Averaging
!!> Average across corners.

subroutine fclaw3dx_clawpatch46_fort_average_corner(mx,my,mz,mbc,meqn, &
    refratio,qcoarse,qfine,areacoarse,areafine, & 
    manifold,icorner_coarse,transform_cptr)
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn,refratio,icorner_coarse, manifold
    INTEGER*8 :: transform_cptr
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold.
    DOUBLE PRECISION :: areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    DOUBLE PRECISION ::   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    INTEGER :: ibc,jbc,mq,r2
    LOGICAL :: is_manifold
    DOUBLE PRECISION :: qf,kf, sum

    !! # This should be refratio*refratio.
    INTEGER :: rr2
    PARAMETER(rr2 = 4)
    INTEGER :: i2(0:rr2-1),j2(0:rr2-1)

    INTEGER :: i1,j1,m, k
    DOUBLE PRECISION :: vf_sum

    r2 = refratio*refratio
    if (r2 .ne. rr2) then
        write(6,*) 'average_corner_ghost (claw2d_utils.f) ', & 
        '  Refratio**2 is not equal to rr2'
        stop
    endif

    is_manifold = manifold .eq. 1

    r2 = refratio*refratio
    !! # Loop over four corner cells on coarse grid
    meqn_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            ibc_loop : do ibc = 1,mbc
                jbc_loop : do jbc = 1,mbc
                    !! # Average fine grid corners onto coarse grid ghost corners
                    if (icorner_coarse .eq. 0) then
                        i1 = 1-ibc
                        j1 = 1-jbc
                    else if (icorner_coarse .eq. 1) then
                        i1 = mx+ibc
                        j1 = 1-jbc
                    else if (icorner_coarse .eq. 2) then
                        i1 = 1-ibc
                        j1 = my+jbc
                    else if (icorner_coarse .eq. 3) then
                        i1 = mx+ibc
                        j1 = my+jbc
                    endif

                    call fclaw3dx_clawpatch_transform_corner_half(i1,j1,i2,j2, transform_cptr)
                    if (is_manifold) then
                        sum = 0
                        vf_sum = 0
                        do m = 0,r2-1
                            qf = qfine(i2(m),j2(m),k,mq)
                            kf = areafine(i2(m),j2(m))
                            sum = sum + kf*qf
                            vf_sum = vf_sum + kf
                        enddo
                        qcoarse(i1,j1,k,mq) = sum/vf_sum
                    else
                        sum = 0
                        do m = 0,r2-1
                            sum = sum + qfine(i2(m),j2(m),k,mq)
                        end do
                    endif
                    qcoarse(i1,j1,k,mq) = sum/dble(r2)
                enddo jbc_loop
            enddo ibc_loop
        enddo k_loop
    end do meqn_loop

end subroutine fclaw3dx_clawpatch46_fort_average_corner


!!> \ingroup  Averaging
!!> Average fine grid siblings to parent coarse grid.
subroutine fclaw3dx_clawpatch46_fort_average2coarse(mx,my,mz,mbc,meqn, & 
           qcoarse,qfine, areacoarse, areafine, igrid,manifold)
    IMPLICIT NONE

    INTEGER :: mx,my,mz,mbc,meqn
    INTEGER :: manifold
    DOUBLE PRECISION :: qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! # these will be empty if we are not on a manifold, and so shouldn't
    !! # be referenced. 
    DOUBLE PRECISION :: areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    DOUBLE PRECISION ::   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    !! # This should be refratio*refratio.
    INTEGER :: rr2
    PARAMETER(rr2 = 4)
    INTEGER :: i2(0:rr2-1),j2(0:rr2-1)


    INTEGER :: i,j,k,mq, ii,jj, i1, j1, r2, m
    INTEGER :: ig, jg, ic_add, jc_add
    INTEGER :: rfactor, refratio, igrid
    LOGICAL :: is_manifold
    DOUBLE PRECISION :: sum, kf, qf, vf_sum

    !! Not sure why there are two different variables here.  Isn't 
    !! rfactor == refratio? 
    rfactor = 2
    refratio = 2

    is_manifold = manifold .eq. 1

    !! # 'iface' is relative to the coarse grid

    r2 = refratio*refratio
    if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ', & 
         '  Refratio**2 is not equal to rr2'
        stop
    endif


    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/rfactor
    jc_add = jg*mx/rfactor

    r2 = refratio*refratio
    meqn_loop : do mq = 1,meqn
        k_loop : do k = 1,mz
            j_loop : do j = 1,my/rfactor
                i_loop : do i = 1,mx/rfactor
                    i1 = i+ic_add
                    j1 = j+jc_add
                    m = 0
                    do jj = 1,refratio
                        do ii = 1,refratio
                           i2(m) = (i-1)*refratio + ii
                           j2(m) = (j-1)*refratio + jj
                           m = m + 1
                        enddo
                    enddo
                    if (is_manifold) then
                        sum = 0
                        do m = 0,r2-1
                            qf = qfine(i2(m),j2(m),k,mq)
                            kf = areafine(i2(m),j2(m))
                            sum = sum + kf*qf
                            vf_sum = vf_sum + kf
                        enddo
                        !! kc = areacoarse(i1,j1,k)
                        qcoarse(i1,j1,k,mq) = sum/vf_sum
                    else
                        sum = 0
                        do m = 0,r2-1
                           qf = qfine(i2(m),j2(m),k,mq)
                           sum = sum + qf
                        enddo
                        qcoarse(i1,j1,k,mq) = sum/r2
                    endif 
                end do i_loop
            enddo j_loop
        enddo k_loop
    enddo meqn_loop
end subroutine  fclaw3dx_clawpatch46_fort_average2coarse


