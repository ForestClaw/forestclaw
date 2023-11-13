!! # check to see if value exceeds threshold

integer function fclaw3d_user_exceeds_th(blockno,& 
                                  qval,qmin,qmax,quad, & 
                                  dx,dy,dz,xc,yc,zc,threshold, &
                                  init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1,-1:1)
    double precision :: dx,dy, dz, xc, yc, zc
    integer :: blockno, init_flag
    logical(kind=4) :: is_ghost
    integer :: refine

    integer*8 cont, fclaw_map_get_context
    integer fclaw2d_map_is_used

    double precision xp, yp, zp

    cont = fclaw_map_get_context()

    if (fclaw2d_map_is_used(cont) .ne. 0) then
        call fclaw_map_3d_c2m(cont,blockno,xc,yc,zc,xp,yp,zp)
    else
        xp = xc
        yp = yc
        zp = zc
    endif


    !! Ghost cells from patch in right half are in region x < 0.5, and so 
    !! would fail tagging criteria
    if (is_ghost) then
        !! Don't use this result
        fclaw3d_user_exceeds_th = -1
        return
    endif

    !! Only refine in right half of domain
    if (yp .lt. 1.0) then
        refine = 1
    else
        refine = 0
    endif

    fclaw3d_user_exceeds_th = refine

end function fclaw3d_user_exceeds_th
