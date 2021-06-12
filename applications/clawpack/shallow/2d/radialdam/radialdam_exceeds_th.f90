!! # check to see if value exceeds threshold

logical function radialdam_exceeds_th(blockno,qval,qmin,qmax,quad, & 
                                      dx,dy,xc,yc,threshhold)
    implicit none
    
    double precision qval,qmin,qmax,threshhold, quad(-1:1,-1:1)
    double precision dx,dy, xc, yc
    integer blockno
    logical refine

    double precision dx2, dy2, dqx,dqy, ds, d, nx, ny
    double precision xp,yp,zp,xpm, xpp, ypm, ypp

    integer*8 cont, get_context
    logical fclaw2d_map_is_used

    cont = get_context()

    IF (fclaw2d_map_is_used(cont)) THEN
        CALL fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
        CALL fclaw2d_map_c2m(cont,blockno,xc-dx,yc,xpm,yp,zp)
        CALL fclaw2d_map_c2m(cont,blockno,xc+dx,yc,xpp,yp,zp)
        CALL fclaw2d_map_c2m(cont,blockno,xc,yc+dy,xp,ypp,zp)
        CALL fclaw2d_map_c2m(cont,blockno,xc,yc-dy,xp,ypm,zp)
        dx2 = xpp - xpm
        dy2 = ypp - ypm
    ELSE
        dx2 = 2*dx
        dy2 = 2*dy
        xp = xc
        yp = yc
    ENDIF

    refine = .false.

    dqx = abs(quad(1,0) - quad(-1,0))/dx2
    dqy = abs(quad(0,1) - quad(0,-1))/dy2
    d = sqrt(xp**2 + yp**2)
    nx = xp/d
    ny = yp/d
    ds = abs(dqx*nx + dqy*ny)
    if (ds .gt. threshhold) then
        refine = .true.
    endif

    radialdam_exceeds_th = refine

end function radialdam_exceeds_th