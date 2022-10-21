double precision function fdisc(blockno,xc,yc,zc)
    use setprob_mod, only : pi, mapping, init_choice, x0, y0, z0, r0, maxelev
    implicit none

    double precision xc,yc, zc
    integer blockno

    double precision xp, yp, zp, rp
    double precision omega0, omega
    double precision z_center, r_sphere_2d

    integer*8 cont, fclaw_map_get_context
    integer fclaw2d_map_is_used

    cont = fclaw_map_get_context()

    if (fclaw2d_map_is_used(cont) .ne. 0) then
        call fclaw3d_map_c2m(cont,blockno,xc,yc,zc,xp,yp,zp)
    else
        xp = xc
        yp = yc
        zp = zc
    endif
    rp = sqrt(xp**2 + yp**2 + zp**2)

    fdisc = 0
    if (init_choice .eq. 1) then
        !! Cylindrical or conical overpressure
        if (mapping .le. 1) then
            !! Initialize data in a cylinder for Cartesian mappings
            fdisc = (xp-x0)**2 + (yp-y0)**2 - r0**2
        else
            !! Cone with central axis p and angle omega0
            !! p = (1,0,0)
            !! d = p(1)*xp + p(2)*yp + p(3)*zp  == xp
            !! To shift to different location, set 'rotate' options.
            omega0 = pi/32.d0
            omega = acos(xp/rp)
            fdisc = omega - omega0            
            !! write(6,*) xp,yp,zp,omega, fdisc
        endif
    else if (init_choice .eq. 2) then
        !! Spherical overpressure 
        if (mapping .le. 1) then
            !! (x0,y0,z0) in Cartesian Coordinates.
            fdisc = (xp-x0)**2 + (yp-y0)**2 + (zp-z0)**2 - r0**2
        else
            !! Radius of 2d sphere
            r_sphere_2d = rp - maxelev*zc

            !! assume sphere is located on z-axis and that z0 in [0,1]
            !! Center of sphere, in physical coordinates.  
            z_center = r_sphere_2d + maxelev*z0

            fdisc = xp**2 + yp**2 + (zp-z_center)**2 - r0**2
        endif
    endif


    return
end function fdisc
