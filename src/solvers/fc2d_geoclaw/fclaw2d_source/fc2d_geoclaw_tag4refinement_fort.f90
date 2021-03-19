SUBROUTINE fc2d_geoclaw_fort_tag4refinement(mx,my,mbc,meqn,maux,xlower,ylower, &
    dx,dy,t,blockno,q,aux,mbathy,level,maxlevel,init_flag,tag_patch)

    USE geoclaw_module, ONLY:dry_tolerance, sea_level
    USE geoclaw_module, ONLY: spherical_distance, coordinate_system

    USE topo_module, ONLY: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    USE topo_module, ONLY: minleveltopo,mtopofiles

    USE topo_module, ONLY: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    USE topo_module, ONLY: minleveldtopo,num_dtopo

    USE qinit_module, ONLY: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
    USE qinit_module, ONLY: min_level_qinit,qinit_type

    USE storm_module, ONLY: storm_specification_type, wind_refine, R_refine, storm_location
    USE storm_module, ONLY: wind_forcing, wind_index, wind_refine

    USE regions_module, ONLY: num_regions, regions, region_type
    USE refinement_module
    IMPLICIT NONE

    INTEGER mx,my, mbc, meqn, tag_patch, init_flag, maux, mbathy
    INTEGER blockno, level, maxlevel
    DOUBLE PRECISION xlower, ylower, dx, dy, t, t0
    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    REAL(kind=8), INTENT(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER i,j, mq,m
    REAL(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    REAL(kind=8) :: speed, eta, ds

    !! Storm specific variables
    REAL(kind=8) :: R_eye(2), wind_speed
    LOGICAL allowflag

    LOGICAL time_interval, space_interval

    real(kind=8) :: xupper, yupper
    INTEGER rmax, tag_patch_regions
    !!logical fc2d_geoclaw_P_intersects_R, region_found
    INTEGER fc2d_geoclaw_refine_using_regions

    tag_patch = 0
    t0 = 0

    xupper = xlower + mx*dx
    yupper = ylower + my*dy

    tag_patch_regions = fc2d_geoclaw_refine_using_regions(level,xlower,ylower,xupper,yupper,t)
    if (tag_patch_regions .ge. 0) then
        !! Tagging based on regions is conclusive
        tag_patch = tag_patch_regions
        return
    endif

    tag_patch = 0  !! Don't refine unless criteria below allows it.

    !! # Refine based only on first variable in system.
    !! Loop over interior points on this grid
    !! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: DO j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy

        x_loop: DO i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

!!            !! Check to see if refinement is forced in any topography file region:
!!            do m=1,mtopofiles
!!                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
!!                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
!!                        y_hi > ylowtopo(m) .AND. y_low < yhitopo(m) ) THEN
!!                        tag_patch = 1
!!                        RETURN
!!                    endif
!!                endif
!!            end do

!!            !! Check to see if refinement is forced in any other region:
!!            DO m=1,num_regions
!!                time_interval = t >= regions(m)%t_low .AND. t <= regions(m)%t_hi
!!                space_interval = x_hi > regions(m)%x_low .AND. x_low < regions(m)%x_hi &
!!                .AND. y_hi > regions(m)%y_low .AND. y_low < regions(m)%y_hi
!!                IF (time_interval .AND. space_interval) THEN
!!                    IF (level < regions(m)%min_level) THEN
!!                        !! Refine to at least to the minimum level
!!                        tag_patch = 1
!!                        RETURN
!!                    ELSEIF (level >= regions(m)%max_level) THEN
!!                        !! Return without refining further ??
!!                    ENDIF
!!                ENDIF
!!            end do

!!            !! This is going to be taken out eventually (5/24/2020)
!!            !! Check if we're in the dtopo region and need to refine:
!!            !! force refinement to level minleveldtopo
!!            do m = 1,num_dtopo
!!                IF (level < minleveldtopo(m).AND. &
!!                    t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
!!                    x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
!!                    y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then
!!                    tag_patch = 1
!!                    RETURN
!!                ENDIF
!!            end do

            !! Check if we're in the region where initial perturbation is
            !! specified and need to force refinement:
            !! This assumes that t0 = 0.d0, should really be t0 but we do
            !! not have access to that parameter in this routine
            ! IF (init_flag .NE. 0) THEN
            IF (qinit_type > 0 .AND. init_flag .NE. 0) THEN
                space_interval = x_hi > x_low_qinit .and. x_low < x_hi_qinit .and. &
                y_hi > y_low_qinit .AND. y_low < y_hi_qinit
                IF (space_interval) THEN
                    IF (level < min_level_qinit) THEN
                        tag_patch = 1
                        RETURN
                    END IF
                endif
            endif

            !! -----------------------------------------------------------------
            !! Refinement not forced, so check if it is allowed and if so,
            !! check if there is a reason to flag this point:
            !!IF (allowflag(x_c,y_c,t,level)) THEN
                if (q(1,i,j) > dry_tolerance) then
                    eta = q(1,i,j) + aux(mbathy,i,j)
                    !! Check wave criteria
                    if (abs(eta - sea_level) > wave_tolerance) then
                        !! Check to see if we are near shore
!!                        IF (q(1,i,j) < deep_depth) THEN
!!                            tag_patch = 1
!!                            return
!!                        ELSE IF (level < max_level_deep) THEN
!!                            tag_patch = 1
!!                            return
!!                        endif
                    endif

                    !! Check speed criteria, note that it might be useful to
                    !! also have a per layer criteria since this is not
                    !! gradient based
                    !! This assumes that mxnest == maxlevel
                    speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                    do m=1,min(size(speed_tolerance),maxlevel)
                        IF (speed > speed_tolerance(m) .AND. level <= m) THEN
                            tag_patch = 1
                            RETURN
                        endif
                    end do
                endif
            !!endif
        end do x_loop
    enddo y_loop

END SUBROUTINE fc2d_geoclaw_fort_tag4refinement


LOGICAL FUNCTION fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,region)
    USE regions_module
    IMPLICIT NONE

    REAL(kind=8) :: xlower, ylower, xupper, yupper, t
    TYPE(region_type) :: region

    if (t < region%t_low .or. t > region%t_hi) then 
        fc2d_geoclaw_P_intersects_R = .false.
        return
    endif

    if (xupper < region%x_low .or. xlower > region%x_hi) then 
        fc2d_geoclaw_P_intersects_R = .false.
        return
    endif

    if (yupper < region%y_low .or. ylower > region%y_hi) then 
        fc2d_geoclaw_P_intersects_R = .false.
        return
    endif

    fc2d_geoclaw_P_intersects_R = .true.
    return

END FUNCTION fc2d_geoclaw_P_intersects_R

INTEGER FUNCTION fc2d_geoclaw_refine_using_regions(level,xlower,ylower,xupper,yupper,t)
    USE regions_module
    IMPLICIT NONE

    REAL(kind=8) :: xlower,ylower,xupper,yupper,t
    integer level

    INTEGER rmax, m, tag_patch, mmax
    LOGICAL region_found, fc2d_geoclaw_P_intersects_R, inregion

    tag_patch = -1  !! use refinement criteria based on wave tolerance, etc

    !! Check to see if refinement is forced by regions R:
    !! If level < R.minlevel for any R : force refinement
    !! If level >= R.maxlevel for all R : don't allow refinement
    !! If none of the above, use usual refinement criteria
    !! If this patch does not intersect any region, use refinement criteria
    rmax = 0
    region_found = .false.
    mmax = 1
    DO m=1,num_regions
        !!inregion = fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,regions(m))
        if (fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,regions(m))) then
            region_found = .true.
            if (level < regions(m)%min_level) then
                !! level < R.minlevel for some R : force refinement
                tag_patch = 1                
            endif
            !! Collect largest max_level
            if (regions(m)%max_level > rmax) then
                mmax = m
                rmax = regions(m)%max_level
            endif
        endif
    end do
    !! level >= R.maxlevel for all R : don't allow refinement
    if (region_found) then
        !! We can use regions as a criteria
        if (level .ge. rmax) then
            !!write(6,'(A,I5,A,I5)') '---------------> level = ',level,'; rmax = ',rmax
            tag_patch = 0
        endif
    endif

    fc2d_geoclaw_refine_using_regions = tag_patch

    return

END FUNCTION fc2d_geoclaw_refine_using_regions
