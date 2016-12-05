SUBROUTINE fc2d_geoclaw_fort_ghostaux(mbc,mx,my,mint,xlow,ylow,dx,dy,maux,aux)
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
  USE geoclaw_module, ONLY: sea_level

  USE storm_module, ONLY: wind_forcing, pressure_forcing
  USE storm_module, ONLY: wind_index, pressure_index, set_storm_fields
  USE storm_module, ONLY: ambient_pressure

  USE friction_module, ONLY: variable_friction, friction_index
  USE friction_module, ONLY: set_ghost_friction_field

  USE topo_module

  IMPLICIT NONE

  !! Arguments
  INTEGER, INTENT(in) :: mbc,mx,my,maux
  REAL(kind=8), INTENT(in) :: xlow,ylow,dx,dy
  REAL(kind=8), INTENT(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  !! For not set value in the center of patch which is not going to be used
  INTEGER :: mint

  !! Locals
  INTEGER :: ii,jj,m,iint,jint
  REAL(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
  CHARACTER(len=*), PARAMETER :: aux_format = "(2i4,4d15.3)"
  INTEGER :: skipcount,iaux,ilo,jlo

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
        
        !!face 0
        aux(2,0:mint,0:my-mint) = 1.d0
        !!face 2
        aux(2,mint+1:mx+1,0:mint) = 1.d0
        !!face 1
        aux(2,mx-mint+1:mx+1,mint+1:my+1) = 1.d0
        !!face 3
        aux(2,0:mx-mint,my-mint+1:my+1) = 1.d0    
        !!aux(3,:,:) = 1.d0
     ENDIF
  ENDIF


  !! If using a variable friction field initialize the coefficients to 0
  IF (variable_friction) THEN
     !!face 0
     aux(friction_index,0:mint,0:my-mint) = 0.d0
     !!face 2
     aux(friction_index,mint+1:mx+1,0:mint) = 0.d0
     !!face 1
     aux(friction_index,mx-mint+1:mx+1,mint+1:my+1) = 0.d0
     !!face 3
     aux(friction_index,0:mx-mint,my-mint+1:my+1) = 0.d0     
  ENDIF

  !! Storm fields if used
  IF (wind_forcing) THEN
     !!face 0
     aux(wind_index,0:mint,0:my-mint) = 0.d0
     aux(wind_index+1,0:mint,0:my-mint) = 0.d0
     !!face 2
     aux(wind_index,mint+1:mx+1,0:mint) = 0.d0
     aux(wind_index+1,mint+1:mx+1,0:mint) = 0.d0
     !!face 1
     aux(wind_index,mx-mint+1:mx+1,mint+1:my+1) = 0.d0
     aux(wind_index+1,mx-mint+1:mx+1,mint+1:my+1) = 0.d0
     !!face 3
     aux(wind_index,0:mx-mint,my-mint+1:my+1) = 0.d0  
     aux(wind_index+1,0:mx-mint,my-mint+1:my+1) = 0.d0    
  ENDIF
  IF (pressure_forcing) THEN
     !!face 0
     aux(pressure_index,0:mint,0:my-mint) = ambient_pressure
     !!face 2
     aux(pressure_index,mint+1:mx+1,0:mint) = ambient_pressure
     !!face 1
     aux(pressure_index,mx-mint+1:mx+1,mint+1:my+1) = ambient_pressure
     !!face 3
     aux(pressure_index,0:mx-mint,my-mint+1:my+1) = ambient_pressure
  ENDIF


  !! Set analytical bathymetry here if requested
  IF (test_topography > 0) THEN
     !!face 0
     FORALL (ii=0:mint,jj=0:my-mint)
        aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
             ylow + (jj - 0.5d0) * dy)
     END FORALL
     !!face 2
     FORALL (ii=mint+1:mx+1,jj=0:mint)
        aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
             ylow + (jj - 0.5d0) * dy)
     END FORALL
     !!face 1
     FORALL (ii=mx-mint+1:mx+1,jj=mint+1:my+1)
        aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
             ylow + (jj - 0.5d0) * dy)
     END FORALL
     !!face 3
     FORALL (ii=0:mx-mint,jj=my-mint+1:my+1)
        aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
             ylow + (jj - 0.5d0) * dy)
     END FORALL
  ENDIF




  !! test:  compute integer indices based off same corner of domain
  !!        to reduce round off discrepancies
  ilo = FLOOR((xlow - xlower + .05d0*dx)/dx)
  jlo = FLOOR((ylow - ylower + .05d0*dy)/dy)

  !! Set bathymetry
  !! Face 0
  DO jj = 0,my-mint
     ym = ylower + (jlo+jj-1.d0) * dy
     yp = ylower + (jlo+jj) * dy
     y  = 0.5d0*(ym+yp)    
     DO ii = 0,mint
        xm = xlower + (ilo+ii-1.d0) * dx
        xp = xlower + (ilo+ii) * dx
        x = 0.5d0*(xm+xp)
!         !!write(*,444)ii,jj,aux(1,ii,jj)
! 444     FORMAT("in setaux ",2i4,e12.5)

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

  !! Face 2
  DO jj = 0,mint
     ym = ylower + (jlo+jj-1.d0) * dy
     yp = ylower + (jlo+jj) * dy
     y  = 0.5d0*(ym+yp)    
     DO ii = mint+1,mx+1
        xm = xlower + (ilo+ii-1.d0) * dx
        xp = xlower + (ilo+ii) * dx
        x = 0.5d0*(xm+xp)
!         !!write(*,444)ii,jj,aux(1,ii,jj)
! 444     FORMAT("in setaux ",2i4,e12.5)

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

  !! Face 1
  DO jj = mint+1,my+1
     ym = ylower + (jlo+jj-1.d0) * dy
     yp = ylower + (jlo+jj) * dy
     y  = 0.5d0*(ym+yp)    
     DO ii = mx-mint+1,mx+1
        xm = xlower + (ilo+ii-1.d0) * dx
        xp = xlower + (ilo+ii) * dx
        x = 0.5d0*(xm+xp)
!         !!write(*,444)ii,jj,aux(1,ii,jj)
! 444     FORMAT("in setaux ",2i4,e12.5)

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

  !! Face 3
  DO jj = my-mint+1,my+1
     ym = ylower + (jlo+jj-1.d0) * dy
     yp = ylower + (jlo+jj) * dy
     y  = 0.5d0*(ym+yp)    
     DO ii = 0,mx-mint
        xm = xlower + (ilo+ii-1.d0) * dx
        xp = xlower + (ilo+ii) * dx
        x = 0.5d0*(xm+xp)
        !!write(*,444)ii,jj,aux(1,ii,jj)
444     FORMAT("in setaux ",2i4,e12.5)

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

  !! Copy topo to ghost cells if outside physical domain
  !! face 0
  DO jj=0,my-mint
     y = ylower + (jlo+jj-.5d0) * dy
     IF ((y < ylower) .OR. (y>yupper)) THEN
        DO ii=0,mint
           x = xlower + (ilo+ii-.5d0) * dx
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO
  !! face 2
  DO jj=0,mint
     y = ylower + (jlo+jj-.5d0) * dy
     IF ((y < ylower) .OR. (y>yupper)) THEN
        DO ii=mint+1,mx+1
           x = xlower + (ilo+ii-.5d0) * dx
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO
  !! face 1
  DO jj=mint+1,my+1
     y = ylower + (jlo+jj-.5d0) * dy
     IF ((y < ylower) .OR. (y>yupper)) THEN
        DO ii=mx-mint+1,mx+1
           x = xlower + (ilo+ii-.5d0) * dx
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO
  !! face 3
  DO jj=my-mint+1,my+1
     y = ylower + (jlo+jj-.5d0) * dy
     IF ((y < ylower) .OR. (y>yupper)) THEN
        DO ii=0,mx-mint
           x = xlower + (ilo+ii-.5d0) * dx
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO

  !! face 0
  DO ii=0,mint
     x =  xlower + (ilo+ii-.5d0) * dx
     IF ((x < xlower) .OR. (x > xupper)) THEN
        DO jj=0,my-mint
           y = ylower + (jlo+jj-.5d0) * dy
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO

  !! face 2
  DO ii=mint+1,mx+1
     x =  xlower + (ilo+ii-.5d0) * dx
     IF ((x < xlower) .OR. (x > xupper)) THEN
        DO jj=0,mint
           y = ylower + (jlo+jj-.5d0) * dy
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO

  !! face 1
  DO ii=mx-mint+1,mx+1
     x =  xlower + (ilo+ii-.5d0) * dx
     IF ((x < xlower) .OR. (x > xupper)) THEN
        DO jj=mint+1,my+1
           y = ylower + (jlo+jj-.5d0) * dy
           iint = ii + MAX(0, CEILING((xlower-x)/dx)) &
                - MAX(0, CEILING((x-xupper)/dx))
           jint = jj + MAX(0, CEILING((ylower-y)/dy)) &
                - MAX(0, CEILING((y-yupper)/dy))
           aux(1,ii,jj) = aux(1,iint,jint)
        ENDDO
     ENDIF
  ENDDO

  !! face 3
  DO ii=0,mx-mint
     x =  xlower + (ilo+ii-.5d0) * dx
     IF ((x < xlower) .OR. (x > xupper)) THEN
        DO jj=my-mint+1,my+1
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
     CALL set_ghost_friction_field(mx,my,mint,mbc,maux,xlow,ylow,dx,dy,aux)
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
END SUBROUTINE fc2d_geoclaw_fort_ghostaux