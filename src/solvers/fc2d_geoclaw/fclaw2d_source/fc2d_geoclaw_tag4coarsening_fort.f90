!  # Template function for setting coarsening criteria.  The
!  # user can copy this file to their directory.  To
!  # indicate that this file should be used, set :
!  #
!  #      fclaw2d_vtable_t vt;
!  #      /* .... */
!  #      vt.fort_tag4coarsening = &tag4coarsening;
!  #      fclaw2d_set_vtable(domain,&vt);
!  #
!  # in virtual tables (typically set in <application>_user.cpp, in a
!  # a routine link '<application>_link_solvers(domain)'
!  #
!  # See also 'tag4refinement.f'

subroutine fc2d_geoclaw_fort_tag4coarsening(blockno,mx,my,mbc,meqn,maux, &
       xlower,ylower,dx,dy,t,q0,q1,q2,q3, &
       aux0,aux1,aux2,aux3,mbathy,level,maxlevel, &
       dry_tolerance_c, wave_tolerance_c, speed_tolerance_entries_c, &
       speed_tolerance_c, tag_patch)

implicit none

INTEGER mx, my, mbc, meqn, maux, tag_patch, blockno, mbathy
INTEGER level, maxlevel, speed_tolerance_entries_c
DOUBLE PRECISION xlower(0:3), ylower(0:3), dx, dy, t, t0
DOUBLE PRECISION wave_tolerance_c, speed_tolerance_c(speed_tolerance_entries_c)
DOUBLE PRECISION dry_tolerance_c

DOUBLE PRECISION q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

REAL(kind=8), INTENT(in) :: aux0(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux1(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux2(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux3(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

call check_patch(mx,my,mbc,meqn,maux,xlower(0),ylower(0), &
                 dx,dy,t,q0,aux0,mbathy,level,maxlevel, &
                 dry_tolerance_c, &
                 wave_tolerance_c,speed_tolerance_entries_c, &
                 speed_tolerance_c, tag_patch)
if (tag_patch == 0) return

call check_patch(mx,my,mbc,meqn,maux,xlower(1),ylower(1), &
                 dx,dy,t,q1,aux1,mbathy,level,maxlevel, &
                 dry_tolerance_c, &
                 wave_tolerance_c,speed_tolerance_entries_c, &
                 speed_tolerance_c, tag_patch)
if (tag_patch == 0) return

call check_patch(mx,my,mbc,meqn,maux,xlower(2),ylower(2), &
                 dx,dy,t,q2,aux2,mbathy,level,maxlevel, &
                 dry_tolerance_c, &
                 wave_tolerance_c,speed_tolerance_entries_c, &
                 speed_tolerance_c, tag_patch)
if (tag_patch == 0) return

call check_patch(mx,my,mbc,meqn,maux,xlower(3),ylower(3), &
                 dx,dy,t,q3,aux3,mbathy,level,maxlevel, &
                 dry_tolerance_c, &
                 wave_tolerance_c,speed_tolerance_entries_c, &
                 speed_tolerance_c, tag_patch)
if (tag_patch == 0) return

return
end subroutine fc2d_geoclaw_fort_tag4coarsening


SUBROUTINE check_patch(mx,my,mbc,meqn,maux,xlower,ylower, &
                       dx,dy,t,q,aux,mbathy,level,maxlevel, &
                       dry_tolerance_c, &
                       wave_tolerance_c,speed_tolerance_entries_c, &
                       speed_tolerance_c,tag_patch)

    USE geoclaw_module, ONLY: sea_level
    USE geoclaw_module, ONLY: spherical_distance, coordinate_system

    USE topo_module, ONLY: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    USE topo_module, ONLY: minleveltopo,mtopofiles

    USE topo_module, ONLY: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    USE topo_module, ONLY: minleveldtopo,num_dtopo

    USE qinit_module, ONLY: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
    USE qinit_module, ONLY: min_level_qinit,qinit_type

    !!USE storm_module, ONLY: storm_specification_type, wind_refine, R_refine, storm_location
    !!USE storm_module, ONLY: wind_forcing, wind_index, wind_refine

    USE regions_module, ONLY: num_regions, regions, region_type
    USE refinement_module

    IMPLICIT NONE

    INTEGER mx, my, mbc, meqn, maux, tag_patch, blockno, mbathy
    INTEGER level, maxlevel, clevel, speed_tolerance_entries_c
    DOUBLE PRECISION xlower, ylower, dx, dy, t, t0
    DOUBLE PRECISION wave_tolerance_c, speed_tolerance_c(speed_tolerance_entries_c)
    DOUBLE PRECISION dry_tolerance_c

    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    REAL(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER i,j,k,m
    REAL(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    REAL(kind=8) :: speed, eta, ds

    !! Storm specific variables
    !!REAL(kind=8) :: R_eye(2), wind_speed

    LOGICAL allowcoarsen

    LOGICAL time_interval, space_interval

    REAL(kind=8) :: xupper, yupper
    integer tag_patch_regions, fc2d_geoclaw_coarsen_using_regions

    tag_patch = 1
    t0 = 1

    !! If coarsened, level will become level - 1
    clevel = level - 1

    xupper = xlower + mx*dx
    yupper = ylower + my*dy

    !! would we refine this patch based on regions criteria? 
    !! If tag_patch_regions == 0 :  level < R.minlevel for some R; must refine --> can't coarsen
    !! If tag_patch_regions == 1 :  level > R.maxlevel for all R; can't refine --> must coarsen
    !! if tag_patch_regions == -1 : Inconclusive (use refinement criteria)
    tag_patch_regions = fc2d_geoclaw_coarsen_using_regions(level,xlower,ylower,xupper,yupper,t)
    if (tag_patch_regions .ge. 0) then
        tag_patch = tag_patch_regions  
        return
    endif

    !! Loop over interior points on this grid
    !! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c   = ylower + (j - 0.5d0) * dy
        y_hi  = ylower + j * dy

        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c   = xlower + (i - 0.5d0) * dx
            x_hi  = xlower + i * dx

            !! Check that the grids is allowed to be coarsened or not
            !!if (allowcoarsen(x_c,y_c,t,clevel)) then
            if (q(1,i,j) > dry_tolerance_c) then
                eta = q(1,i,j) + aux(mbathy,i,j)
                if ( abs(eta - sea_level) > wave_tolerance_c) then
                    !! This patch fails the coarsening test
                    tag_patch = 0
                else
                    tag_patch = 1
                    return
                endif
                !! Ignore the speed criteria first
                ! speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                ! if ( abs(eta - sea_level) < wave_tolerance_c &
                !     .and. speed < speed_tolerance_c(level)) then
                !     tag_patch = 1
                ! else
                !     tag_patch = 0
                !     return
                ! endif
            else
                !! we can always coarsen if we are on dry land
                tag_patch = 1
            endif
!!     else
!!        tag_patch = 0
!!        return
     !!endif

        end do x_loop
    end do y_loop

END SUBROUTINE check_patch


!! This looks almost like the refine version
!! In the refine version, we don't allow refinement if level >= rmax
!! In the coarsen version we allow coarsening only if level > rmax
integer function fc2d_geoclaw_coarsen_using_regions(level,xlower,ylower,xupper,yupper,t)
    USE regions_module
    IMPLICIT NONE

    REAL(kind=8) :: xlower,ylower,xupper,yupper,t
    integer level

    INTEGER rmax, m, tag_patch
    LOGICAL region_found, fc2d_geoclaw_P_intersects_R

    tag_patch = -1  !! inconclusive

    !! Check to see if refinement is forced by regions :
    !! If level < R.minlevel for any R : force refinement
    !! If level >= R.maxlevel for all R : don't allow refinement
    rmax = 0
    region_found = .false.
    DO m=1,num_regions
        if (fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,regions(m))) then
            region_found = .true.
            if (level < regions(m)%min_level) then
                !! level < R.minlevel for any R : force refinement (don't coarsen)
                tag_patch = 0
                return
            endif
            !! Collect largest max_level
            if (regions(m)%max_level > rmax) then
                rmax = regions(m)%max_level
            endif
        endif
    end do
    !! level >= R.maxlevel for all R : don't allow refinement
    if (region_found) then
        !! We can use regions as a criteria
        if (level .gt. rmax) then
            !! Force coarsening
            tag_patch = 1
        endif
    endif

    fc2d_geoclaw_coarsen_using_regions = tag_patch

    return
END FUNCTION fc2d_geoclaw_coarsen_using_regions


logical function allowcoarsen(x,y,t,level)

use amr_module, only: t0
use geoclaw_module
use regions_module
use refinement_module
use topo_module
use qinit_module

implicit none

! Function arguments
real(kind=8), intent(in) :: x,y,t
integer, intent(in) :: level

! Locals
integer :: m

allowcoarsen = .true.

!   following commented by dlg on 10/9/08.
!   my understanding of maxleveldeep might be differnet
!   still shouldn't be allowed if maxlevel allowed in a region is less
!   than maxleveldeep
!   might want to allow high levels of refinement in some deep regions
!   but not others.
!
!   if (level .lt. maxleveldeep) then
!      # allow refinement to next level in deep water
!       allowflag = .true.
!       go to 900  !# no need to check anything else
!       endif

!! t should be t + dt_coarsened
! Allow flagging everywhere if using test bathymetry
! if (test_topography > 1) then
!     allowcoarsen = .true.
!     return
! endif

do m=1,mtopofiles
  if (x > xlowtopo(m) .and. x < xhitopo(m) .and. &
      y > ylowtopo(m) .and. y < yhitopo(m) .and. &
      t >= tlowtopo(m) .and. t < thitopo(m)) then
      if (level < minleveltopo(m)) then
        allowcoarsen = .false.
        return
      endif
  endif
enddo

do m=1,num_regions
  if (level < regions(m)%min_level) then
    if (x > regions(m)%x_low .and. x <  regions(m)%x_hi.and. &
        y > regions(m)%y_low .and. y <  regions(m)%y_hi.and. &
        t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then
        allowcoarsen = .false.
        return
      endif
  endif
enddo

do m=1,num_dtopo
    if (x >  xlowdtopo(m) .and. x < xhidtopo(m).and. &
        y >  ylowdtopo(m) .and. y < yhidtopo(m).and. &
        t >= t0dtopo(m)   .and. t <= tfdtopo(m)) then
        if (level < minleveldtopo(m)) then
            allowcoarsen = .false.
            return
        endif
    endif
enddo

end function allowcoarsen
