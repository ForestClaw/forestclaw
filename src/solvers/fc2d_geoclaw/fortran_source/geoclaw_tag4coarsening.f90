! c     # Template function for setting coarsening criteria.  The
! c     # user can copy this file to their directory.  To
! c     # indicate that this file should be used, set :
! c     #
! c     #      fclaw2d_vtable_t vt;
! c     #      /* .... */
! c     #      vt.fort_tag4coarsening = &tag4coarsening;
! c     #      fclaw2d_set_vtable(domain,&vt);
! c     #
! c     # in virtual tables (typically set in <application>_user.cpp, in a
! c     # a routine link '<application>_link_solvers(domain)'
! c     #
! c     # See also 'tag4refinement.f'

subroutine geoclaw_tag4coarsening(mx,my,mbc,meqn,maux, &
       xlower,ylower,dx,dy,blockno,q0,q1,q2,q3, &
       aux0,aux1,aux2,aux3, level, maxlevel, tag_patch)

USE geoclaw_module, ONLY:dry_tolerance, sea_level
USE geoclaw_module, ONLY: spherical_distance, coordinate_system

USE topo_module, ONLY: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
USE topo_module, ONLY: minleveltopo,mtopofiles

USE topo_module, ONLY: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
USE topo_module, ONLY: minleveldtopo,num_dtopo

USE qinit_module, ONLY: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
USE qinit_module, ONLY: min_level_qinit,qinit_type

USE storm_module, ONLY: storm_type, wind_refine, R_refine, storm_location
USE storm_module, ONLY: wind_forcing, wind_index, wind_refine

USE regions_module, ONLY: num_regions, regions, region_type
USE refinement_module
implicit none

INTEGER mx,my, mbc, meqn, maux, tag_patch
INTEGER blockno, level, maxlevel
DOUBLE PRECISION xlower(0:3), ylower(0:3), dx, dy, t, t0

DOUBLE PRECISION q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
DOUBLE PRECISION q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

REAL(kind=8), INTENT(in) :: aux0(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux1(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux2(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8), INTENT(in) :: aux3(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
REAL(kind=8) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

INTEGER i,j,k,m
REAL(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
REAL(kind=8) :: speed, eta, ds

!! Storm specific variables
REAL(kind=8) :: R_eye(2), wind_speed
LOGICAL allowcoarsen

LOGICAL time_interval, space_interval

tag_patch = 1
t0 = 0

level = level - 1
!! # Refine based only on first variable in system.
!! Loop over interior points on this grid
!! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
q_loop: do k=1,4
   select case(k)
      case (0)
         q = q0
         aux = aux0
      case (1)
         q = q1
         aux = aux1
      case (2)
         q = q2
         aux = aux2
      case (3)
         q = q3
         aux = aux3
   end select

   y_loop: do j=1-mbc,my+mbc
     y_low = ylower(k) + (j - 1) * dy
     y_c = ylower(k) + (j - 0.5d0) * dy
     y_hi = ylower(k) + j * dy

     x_loop: do i = 1-mbc,mx+mbc
        x_low = xlower(k) + (i - 1) * dx
        x_c = xlower(k) + (i - 0.5d0) * dx
        x_hi = xlower(k) + i * dx 

        !! (Ignore the storm refinement first.)

        !! Here, t should be t + dt_coarsened (Not solved yet)
        !! Check to see if refinement is forced in any topography file region:
        do m=1,mtopofiles
           if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
              if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                    y_hi > ylowtopo(m) .AND. y_low < yhitopo(m) ) THEN
                 tag_patch = 0
                 RETURN
              endif
           endif
        enddo

        !! Check to see if refinement is forced in any other region:
        DO m=1,num_regions
           time_interval = t >= regions(m)%t_low .AND. t <= regions(m)%t_hi
           space_interval = x_hi > regions(m)%x_low .AND. x_low < regions(m)%x_hi &
                .AND. y_hi > regions(m)%y_low .AND. y_low < regions(m)%y_hi
           IF (time_interval .AND. space_interval) THEN
              IF (level < regions(m)%min_level) THEN
                 !! Refine to at least to the minimum level
                 tag_patch = 0
                 RETURN
              ENDIF
           ENDIF
        enddo

        !! Check if we're in the dtopo region and need to refine:
        !! force refinement to level minleveldtopo
        do m = 1,num_dtopo
           IF (level < minleveldtopo(m).AND. &
                t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
                x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
                y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then
              tag_patch = 0
              RETURN
           ENDIF
        enddo

        !! -----------------------------------------------------------------
        !! Refinement not forced, so check if it is allowed and if so,
        !! check if there is a reason to flag this point:
        if (allowcoarsen(x_c,y_c,t,level)) then
           if (q(1,i,j) > dry_tolerance) then
              eta = q(1,i,j) + aux(1,i,j)
              speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)

              if ( abs(eta - sea_level) < wave_tolerance &
                  .and. speed < speed_tolerance(level)) then
                  tag_patch = 1
              else 
                  tag_patch = 0
                  return
              endif
              !! Check wave criteria
              ! if (abs(eta - sea_level) > wave_tolerance) then
              !    !! Check to see if we are near shore
              !    if (q(1,i,j) < deep_depth) then
              !       tag_patch = 0
              !       return
              !       !! Check if we are allowed to flag in deep water
              !       !! anyway
              !    else if (level < max_level_deep) then
              !       tag_patch = 0
              !       return
              !    endif
              ! endif

              !! Check speed criteria, note that it might be useful to
              !! also have a per layer criteria since this is not
              !! gradient based
              !! This assumes that mxnest == maxlevel
              ! speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
              ! do m=1,min(size(speed_tolerance),maxlevel)
              !    IF (speed > speed_tolerance(m) .AND. level <= m) THEN
              !       tag_patch = 0
              !       RETURN
              !    endif
              ! enddo
           endif
        endif

     enddo x_loop
   enddo y_loop
enddo q_loop

end

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

    allowcoarsen = .false.

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
    if (test_topography > 1) then
        allowcoarsen = .true.
        return
    endif
    do m=1,mtopofiles
        if (level > minleveltopo(m)) then
            if (x > xlowtopo(m) .and. x < xhitopo(m) .and. &
                y > ylowtopo(m) .and. y < yhitopo(m) .and. &
                t >= tlowtopo(m) .and. t < thitopo(m)) then

                allowcoarsen = .true.
                return
            endif
        endif
    enddo
    do m=1,num_regions
        if (level > regions(m)%min_level) then
            if (x > regions(m)%x_low .and. x <  regions(m)%x_hi.and. &
                y > regions(m)%y_low .and. y <  regions(m)%y_hi.and. &
                t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then

                allowcoarsen = .true.
                return
            endif
        endif
    enddo

    do m=1,num_dtopo
        if (x >  xlowdtopo(m) .and. x < xhidtopo(m).and. &
            y >  ylowdtopo(m) .and. y < yhidtopo(m).and. &
            t >= t0dtopo(m)   .and. t <= tfdtopo(m)) then

            if (level > minleveldtopo(m)) then
                allowcoarsen = .true.
                return
            endif
        endif
    enddo

end function allowcoarsen
