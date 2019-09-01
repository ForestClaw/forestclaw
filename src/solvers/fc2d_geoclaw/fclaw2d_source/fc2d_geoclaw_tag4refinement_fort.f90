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

  tag_patch = 0
  t0 = 0

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

        !! The following conditions are only checked in the horizontal and
        !! override the allowflag routine

        !! ************* Storm Based Refinement ****************
        !! Check to see if we are some specified distance from the eye of
        !! the storm and refine if we are
        if (storm_specification_type > 0) then
           R_eye = storm_location(t)
           do m=1,size(R_refine,1)
              if (coordinate_system == 2) then
                 ds = spherical_distance(x_c, y_c, R_eye(1), R_eye(2))
              else
                 ds = sqrt((x_c - R_eye(1))**2 + (y_c - R_eye(2))**2)
              end if

              IF ( ds < R_refine(m) .AND. level <= m ) THEN
                 tag_patch = 1
                 return
              endif
           enddo

           !! Refine based on wind speed
           if (wind_forcing) then
              wind_speed = sqrt(aux(wind_index,i,j)**2 + aux(wind_index+1,i,j)**2)
              do m=1,size(wind_refine,1)
                 IF ((wind_speed > wind_refine(m)) .AND. (level <= m)) THEN
                    tag_patch = 1
                    RETURN
                 endif
              enddo
           endif
        endif
        !! *****************************************************

        !! Check to see if refinement is forced in any topography file region:
        do m=1,mtopofiles
           if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
              if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                    y_hi > ylowtopo(m) .AND. y_low < yhitopo(m) ) THEN
                 tag_patch = 1
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
                 tag_patch = 1
                 RETURN
              ELSEIF (level >= regions(m)%max_level) THEN
                 !! Return without refining further ??
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
              tag_patch = 1
              RETURN
           ENDIF
        enddo

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
        IF (allowflag(x_c,y_c,t,level)) THEN
           if (q(1,i,j) > dry_tolerance) then
              eta = q(1,i,j) + aux(mbathy,i,j)
              !! Check wave criteria
              if (abs(eta - sea_level) > wave_tolerance) then
                 !! Check to see if we are near shore
                 IF (q(1,i,j) < deep_depth) THEN
                    tag_patch = 1
                    ! write (*,*) 'near shore: x_c, y_c, t, level', x_c,y_c,t,level
                    ! write (*,*) 'perturbation, wave tolerance', abs(eta - sea_level), wave_tolerance
                    return
                    !! Check if we are allowed to flag in deep water
                    !! anyway
                 ELSE IF (level < max_level_deep) THEN
                    tag_patch = 1
                    ! write (*,*) 'not in deep water: x_c, y_c, t, level', x_c,y_c,t,level
                    return
                 endif
              endif

              !! Check speed criteria, note that it might be useful to
              !! also have a per layer criteria since this is not
              !! gradient based
              !! This assumes that mxnest == maxlevel
              speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
              do m=1,min(size(speed_tolerance),maxlevel)
                 IF (speed > speed_tolerance(m) .AND. level <= m) THEN
                    tag_patch = 1
                    ! write (*,*) 'speed: x_c, y_c, t, sea_levelel', x_c,y_c,t,level
                    RETURN
                 endif
              enddo
           endif
        endif

     enddo x_loop
  enddo y_loop

END SUBROUTINE fc2d_geoclaw_fort_tag4refinement
