SUBROUTINE fc2d_geoclaw_fort_tag4coarsening(blockno,mx,my,mbc,meqn,maux, &
       xlower,ylower,dx,dy,t,q0,q1,q2,q3, &
       aux0,aux1,aux2,aux3,level,maxlevel, init_flag, tag_patch)

    IMPLICIT NONE

    INTEGER :: mx, my, mbc, meqn, maux, tag_patch, blockno
    INTEGER :: level, maxlevel, init_flag
    DOUBLE PRECISION :: xlower(0:3), ylower(0:3), dx, dy, t

    DOUBLE PRECISION :: q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    DOUBLE PRECISION :: aux0(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: aux1(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: aux2(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: aux3(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

!!    !! # Don't coarsen when initializing the mesh
!!    if (initflag .ne. 0) then
!!        tag_patch = 0
!!        return
!!    endif

    !! Assume we will coarsen, unless one of the four patches fails
    !! the coarsening test.
    tag_patch = 1

    call fc2d_geoclaw_test_coarsen(blockno, mx,my,mbc,meqn,maux, & 
                xlower(0), ylower(0), & 
                dx,dy,t,q0,aux0,level,maxlevel, init_flag, tag_patch)
    if (tag_patch == 0) return  !!  Can't coarsen

    call fc2d_geoclaw_test_coarsen(blockno, mx,my,mbc,meqn,maux, & 
                xlower(1), ylower(1), & 
                dx,dy,t,q1,aux1,level,maxlevel, init_flag, tag_patch)
    if (tag_patch == 0) return  !!  Can't coarsen  

    call fc2d_geoclaw_test_coarsen(blockno, mx,my,mbc,meqn,maux, &
                xlower(2), ylower(2), & 
                 dx,dy,t,q2,aux2,level,maxlevel, init_flag, tag_patch)
    if (tag_patch == 0) return  !!  Can't coarsen

    call fc2d_geoclaw_test_coarsen(blockno, mx,my,mbc,meqn,maux,&
                xlower(3), ylower(3), & 
                 dx,dy,t,q3,aux3,level,maxlevel, init_flag, tag_patch)
    if (tag_patch == 0) return  !!  Can't coarsen

    return
end subroutine fc2d_geoclaw_fort_tag4coarsening


SUBROUTINE fc2d_geoclaw_test_coarsen(blockno, mx,my,mbc,meqn,maux,xlower,ylower, &
                       dx,dy,t,q,aux,level,maxlevel, init_flag, tag_patch)

    IMPLICIT NONE

    INTEGER :: mx, my, mbc, meqn, maux, tag_patch, blockno
    INTEGER :: level, maxlevel,init_flag
    DOUBLE PRECISION :: xlower, ylower, dx, dy, t

    DOUBLE PRECISION :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER :: i,j,m
    DOUBLE PRECISION :: xc,yc, xupper, yupper, qvec(meqn), auxvec(maux)

    !INTEGER :: tag_patch_regions!, fc2d_geoclaw_coarsen_using_regions
    LOGICAL :: is_coarsening

    INTEGER :: fc2d_geoclaw_flag2refine, flag_patch


    !! Check if regions would allow coarsening
    xupper = xlower + mx*dx
    yupper = ylower + my*dy

    !! If tag_patch_regions == 0 :  level < R.minlevel for some R; must refine --> can't coarsen
    !! If tag_patch_regions == 1 :  level > R.maxlevel for all R; can't refine --> must coarsen
    !! if tag_patch_regions == -1 : Inconclusive (use refinement criteria)

!!    tag_patch_regions = fc2d_geoclaw_coarsen_using_regions(level,xlower,ylower,xupper,yupper,t)
!!    if (tag_patch_regions .ge. 0) then
!!        !! We have to coarsen
!!        tag_patch = tag_patch_regions  
!!        return
!!    endif

    tag_patch = 1    !! Allow coarsening if nothing below prevents it

    is_coarsening = .true.
    do j = 1,my
        yc   = ylower + (j - 0.5d0) * dy
        do i = 1,mx
            xc   = xlower + (i - 0.5d0) * dx

            do m = 1,meqn
                qvec(m) = q(m,i,j)
            end do
            do m = 1,maux
                auxvec(m) = aux(m,i,j)
            enddo

            flag_patch = fc2d_geoclaw_flag2refine( & 
                    blockno, meqn, maux, qvec, auxvec, dx,dy,xc,yc,t,level, & 
                    maxlevel, init_flag, is_coarsening)

!!          # flag_patch : 
!!          # -1 : Not conclusive (possibly ghost cell) (do not tag for coarsening)
!!          # 0  : Does not exceed threshold (tag for coarsening)      
!!          # 1  : Exceeds coarsening threshold (do not tag for coarsening)

            if (flag_patch .gt. 0) then
                tag_patch = 0
                return
            endif

        end do
    end do

END SUBROUTINE fc2d_geoclaw_test_coarsen


!! This looks almost like the refine version
!! In the refine version, we don't allow refinement if level >= rmax
!! In the coarsen version we allow coarsening only if level > rmax

INTEGER FUNCTION fc2d_geoclaw_coarsen_using_regions(level,xlower,ylower,xupper,yupper,t)
    USE regions_module
    IMPLICIT NONE

    double precision :: xlower,ylower,xupper,yupper,t
    integer :: level

    INTEGER :: rmax, m, tag_patch, mmax
    LOGICAL :: region_found, fc2d_geoclaw_P_intersects_R

    tag_patch = -1  !! inconclusive

    !! Check to see if refinement is forced by regions :
    !! If level < R.minlevel for any R : force refinement
    !! If level >= R.maxlevel for all R : don't allow refinement
    rmax = 0
    region_found = .false.
    DO m = 1,num_regions
        if (fc2d_geoclaw_P_intersects_R(xlower,ylower,xupper,yupper,t,regions(m))) then
            region_found = .true.
            if (level < regions(m)%min_level) then
                !! level < R.minlevel for any R : force refinement (don't coarsen)
                tag_patch = 0
                return
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
        if (level .gt. rmax) then
            !! Force coarsening
            tag_patch = 1
        endif
    endif

    fc2d_geoclaw_coarsen_using_regions = tag_patch

    RETURN
END FUNCTION fc2d_geoclaw_coarsen_using_regions


LOGICAL FUNCTION allowcoarsen(x,y,t,level)

!use amr_module, only: t0 
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
!!      if (level < minleveltopo(m)) then
!!        allowcoarsen = .false.
!!        return
!!      endif
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
!!        if (level < minleveldtopo(m)) then
!!            allowcoarsen = .false.
!!            return
!!        endif
    endif
enddo

end function allowcoarsen
