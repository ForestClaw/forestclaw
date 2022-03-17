!
SUBROUTINE fc2d_geoclaw_setaux(mbc,mx,my,xlow,ylow,dx,dy,maux,aux, & 
                               is_ghost_in,nghost,mint)
  !!     ============================================
  !!
  !!     # set auxiliary arrays
  !!
  !!     aux(1,i,j) = Z(x,y) topography (negative below sea level for topoymetry)
  !!
  !!     If coordinate_system=2 then lat-lon coordinates on the sphere and
  !!        aux(2,i,j) = area ratio (capacity function -- set mcapa = 2)
  !!        aux(3,i,j) = length ratio for edge
  !!
  !!

  USE amr_module, ONLY: mcapa, xupper, yupper, xlower, ylower, NEEDS_TO_BE_SET

  USE geoclaw_module, ONLY: coordinate_system, earth_radius, deg2rad
  USE geoclaw_module, ONLY: sea_level, ambient_pressure

  USE storm_module, ONLY: wind_forcing, pressure_forcing
  USE storm_module, ONLY: wind_index, pressure_index, set_storm_fields

  USE friction_module, ONLY: variable_friction, friction_index
  USE friction_module, ONLY: set_friction_field

  USE topo_module

  IMPLICIT NONE

  !! Arguments
  INTEGER, INTENT(in) :: mbc,mx,my,maux
  REAL(kind=8), INTENT(in) :: xlow,ylow,dx,dy
  REAL(kind=8), INTENT(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  integer, intent(in) :: nghost, mint
  INTEGER, INTENT(in) :: is_ghost_in
  LOGICAL*1 :: is_ghost

  !! Locals
  INTEGER :: ii,jj,m, iint,jint
  REAL(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
  CHARACTER(len=*), PARAMETER :: aux_format = "(2i4,4d15.3)"
  INTEGER :: skipcount,iaux,ilo,jlo
  LOGICAL ghost_invalid

  is_ghost = is_ghost_in .ne. 0


  !! Lat-Long coordinate system in use, check input variables
  IF (coordinate_system == 2) THEN
     IF (mcapa /= 2 .OR. maux < 3) THEN
        PRINT *,'ERROR in setaux:  for coordinate_system==2'
        PRINT *,'     need mcapa == 2 and maux >= 3'
        PRINT *,'     have mcapa = ',mcapa,'  maux = ',maux
        STOP
     ENDIF
  ENDIF

  !! Check below is new in 5.2.1 -- need to rethink for storm surge
  !! and other applications where other aux arrays might be used?
  IF (coordinate_system == 1) THEN
     IF (mcapa > 0) THEN
        PRINT *,'ERROR in setaux:  for coordinate_system==1'
        PRINT *,'     need mcapa == 0 and maux == 1'
        PRINT *,'     have mcapa = ',mcapa,'  maux = ',maux
        STOP
     ELSE IF (maux > 1) THEN
        !! Should not need to set aux(2,:,:) in this case, but
        !! for some reason it bombs, e.g. in bowl-radial if maux>1.
        do jj = 1-mbc,my+mbc
          do ii = 1-mbc,mx+mbc
            if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
              cycle
            endif
            aux(2,ii,jj) = 1.d0
          enddo
        enddo
     ENDIF
  ENDIF


  !! If using a variable friction field initialize the coefficients to 0
  IF (variable_friction) THEN
    do jj = 1-mbc,my+mbc
      do ii = 1-mbc,mx+mbc
        if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
          cycle
        endif
        aux(friction_index,ii,jj) = 0.d0
      enddo
    enddo
  ENDIF

    !! Storm fields if used
  IF (wind_forcing) THEN
    do jj = 1-mbc,my+mbc
      do ii = 1-mbc,mx+mbc
        if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
          cycle
        endif
        aux(wind_index,ii,jj) = 0.d0
        aux(wind_index + 1,ii,jj) = 0.d0
      enddo
    enddo
  ENDIF

  IF (pressure_forcing) THEN
    do jj = 1-mbc,my+mbc
      do ii = 1-mbc,mx+mbc
        if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
          cycle
        endif
        aux(pressure_index,ii,jj) = ambient_pressure
      enddo
    enddo
  ENDIF

  !! Set analytical bathymetry here if requested
  IF (test_topography > 0) THEN
    do jj = 1-mbc,my+mbc
      do ii = 1-mbc,mx+mbc
        if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
          cycle
        endif
!!        aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
!!            ylow + (jj - 0.5d0) * dy)
          aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx)
      enddo
    enddo
  ENDIF

  !! test:  compute integer indices based off same corner of domain
  !!        to reduce round off discrepancies
  ilo = FLOOR((xlow - xlower + .05d0*dx)/dx)
  jlo = FLOOR((ylow - ylower + .05d0*dy)/dy)

  !! Set bathymetry
  skipcount = 0
  DO jj=1-mbc,my+mbc
     !!ym = ylow + (jj - 1.d0) * dy
     !!y = ylow + (jj - 0.5d0) * dy
     !!yp = ylow + real(jj,kind=8) * dy

     ym = ylower + (jlo+jj-1.d0) * dy
     yp = ylower + (jlo+jj) * dy
     y = 0.5d0*(ym+yp)


     DO ii=1-mbc,mx+mbc
        !!xm = xlow + (ii - 1.d0) * dx
        !!x  = xlow + (ii - 0.5d0) * dx
        !!xp = xlow + real(ii,kind=8) * dx

        xm = xlower + (ilo+ii-1.d0) * dx
        xp = xlower + (ilo+ii) * dx
        x = 0.5d0*(xm+xp)



        !!write(*,444)ii,jj,aux(1,ii,jj)
444     FORMAT("in setaux ",2i4,e12.5)

        IF (is_ghost .AND. ghost_invalid(ii,jj,mx,my,nghost,mint)) THEN
           CYCLE
        ENDIF


        !! Set lat-long cell info
        IF (coordinate_system == 2) THEN
           aux(2,ii,jj) = deg2rad * earth_radius**2 * (SIN(yp * deg2rad) - SIN(ym * deg2rad)) / dy
           aux(3,ii,jj) = ym * deg2rad
        ENDIF

        !! skip setting aux(1,ii,jj) in ghost cell if outside physical domain
        !! since topo files may not cover ghost cell, and values
        !! should be extrapolated, which is done in next set of loops.
        IF ((y>yupper) .OR. (y<ylower) .OR. &
             (x>xupper) .OR. (x<xlower)) CYCLE

#if 0
        !! ### parameter NEEDS_TO_BE_SET initialized in amr_module.f90
        !! ### saves time by otherwise copying instead of reinitializing
        !! FORESTCLAW change
        IF (aux(1,ii,jj) .NE. NEEDS_TO_BE_SET) THEN
           skipcount = skipcount + 1
           CYCLE  ! new system copies bathy where possible
        ENDIF
#endif



        !! Use input topography files if available
        IF (mtopofiles > 0 .AND. test_topography == 0) THEN
           topo_integral = 0.d0
           CALL cellgridintegrate(topo_integral,xm,x,xp,ym,y,yp, &
                xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                mxtopo,mytopo,mtopo,i0topo,mtopoorder, &
                mtopofiles,mtoposize,topowork)

           IF (coordinate_system == 2) THEN
              aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
           ELSE
              aux(1,ii,jj) = topo_integral / (dx * dy)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !!write(*,*)" skipcount = ",skipcount

  !! Copy topo to ghost cells if outside physical domain

  DO jj=1-mbc,my+mbc
     y = ylower + (jlo+jj-.5d0) * dy
     IF ((y < ylower) .OR. (y>yupper)) THEN
        DO ii=1-mbc,mx+mbc
           if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
             cycle
           endif
           x = xlower + (ilo+ii-.5d0) * dx
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO


  DO ii=1-mbc,mx+mbc
     x =  xlower + (ilo+ii-.5d0) * dx
     IF ((x < xlower) .OR. (x > xupper)) THEN
        DO jj=1-mbc,my+mbc
           if (is_ghost .and. ghost_invalid(ii,jj,mx,my,nghost,mint)) then
             cycle
           endif
           y = ylower + (jlo+jj-.5d0) * dy
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO


  !! Set friction coefficient based on a set of depth levels
  IF (friction_index > 0) THEN
     CALL set_friction_field(mx,my,mbc,maux,xlow,ylow,dx,dy,aux &
                             ,is_ghost,nghost,mint)
  ENDIF


  !! Output for debugging to fort.23
  IF (.FALSE.) THEN
     PRINT *,'Writing out aux arrays'
     PRINT *,' '
     WRITE(23,230)  mbc,mx,my,dx,dy,xlow,ylow
230  FORMAT('==> mbc, mx, my:  ',3i5,'  dx, dy:',2f10.6, &
          '  xlow,ylow:', 2f10.6)
     DO jj=1-mbc,my+mbc
        DO ii=1-mbc,mx+mbc
           x = xlow + (ii-0.5d0)*dx
           y = ylow + (jj-0.5d0)*dy
           IF ((x>223) .AND. (x<232) .AND. (y<37)) &
                WRITE(23,231) ii,jj,x,y,(aux(m,ii,jj),m=1,maux)
231        FORMAT(2i4,2f10.3,3e20.10)
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE fc2d_geoclaw_setaux

!! Change this to ghost_valid
LOGICAL FUNCTION ghost_invalid(i,j,mx,my,nghost,mint)
  implicit none
  integer, intent(in) :: i,j,nghost,mint,mx,my
  logical :: inner, outer

  inner = (i .gt. mint .and. i .lt. mx-mint+1) .and. &
          (j .gt. mint .and. j .lt. my-mint+1)

  outer = (i .lt. 1-nghost) .or. (i .gt. mx+nghost) .or. &
          (j .lt. 1-nghost) .or. (j .gt. my+nghost)

  ghost_invalid = (inner .or. outer)
end function ghost_invalid
