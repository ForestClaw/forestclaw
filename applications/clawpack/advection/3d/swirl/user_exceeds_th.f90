!! # check to see if value exceeds threshold

integer function fclaw3dx_user_exceeds_th(blockno,& 
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

    !! Ghost cells from patch in right half are in region x < 0.5, and so 
    !! would fail tagging criteria
    if (is_ghost) then
        !! Don't use this result
        fclaw3dx_user_exceeds_th = -1
        return
    endif

    !! Only refine in right half of domain
    if (xc .gt. 0.5) then
        refine = 1
    else
        refine = 0
    endif

    fclaw3dx_user_exceeds_th = refine

end function fclaw3dx_user_exceeds_th
