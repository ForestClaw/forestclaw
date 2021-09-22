!! Determine tagging based on regions
SUBROUTINE fc2d_geoclaw_test_regions(level,xlower,ylower,xupper,yupper, & 
                                     t,refine, tag_patch)
    USE regions_module
    IMPLICIT NONE

    DOUBLE PRECISION :: xlower,ylower,xupper,yupper,t
    integer :: level, refine, tag_patch

    INTEGER :: m, min_level, max_level
    LOGICAL :: region_found, fc2d_geoclaw_P_intersects_R

    LOGICAL :: inregion(num_regions)

    tag_patch = -1  !!  Inconclusive for now.

    !! Get list of regions that this patch intersects
    !! If we are coarsening, the "patch" dimensions are the dimensions of the 
    !! quadrant occupied by parent quadrant, i.e. the coarsened patch.  But 'level'
    !! is the level of the four siblings.
    region_found = .false.
    DO m = 1,num_regions
        inregion(m) = fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,regions(m))
        if (inregion(m)) then
            region_found = .true.
        endif
    end do
    if (.not. region_found) then
        !! Refinement criteria not be based on regions
        tag_patch = -1
        return
    endif

    !! Find minimum and maximum levels for regions intersected by this patch
    min_level = 100    !! larger than any possible number of levels
    max_level = 0
    do m = 1,num_regions
        if (inregion(m)) then            
            min_level = min(min_level,regions(m)%min_level)
            max_level = max(max_level,regions(m)%max_level)
        endif
    end do

    !! Determine if we are allowed to refine or coarsen, based on regions above.
    if (refine .ne. 0) then
        !! We are tagging for refinement
        if (level .lt. min_level) then
            !! We have to refine
            tag_patch = 1
        elseif (level .ge. max_level) then
            !! We are not allowed to refine
            tag_patch = 0
        endif
    else
        !! We are tagging for coarsening
        if (level .le. min_level) then
            !! We can't coarsen
            tag_patch = 0
        elseif (level .gt. max_level) then
            !! We have to coarsen
            tag_patch = 1
        endif
    endif

    return

END SUBROUTINE fc2d_geoclaw_test_regions


!! This the old style region tagging and does not yet check
!! the rectangular regions from 'flagregions.f90' in geoclaw.
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