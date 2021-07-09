!! # check to see if value exceeds threshold

integer function fclaw2d_clawpatch_user_exceeds_th(blockno,& 
                                  qval,qmin,qmax,quad, & 
                                  dx,dy,xc,yc,threshold, &
                                  init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag
    logical(kind=4) :: is_ghost
    integer :: refine

    refine = 0
    if (qval .gt. threshold) then
        refine = 1
    endif

    fclaw2d_clawpatch_user_exceeds_th = refine

end function fclaw2d_clawpatch_user_exceeds_th
