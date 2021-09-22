!! # check to see if value exceeds threshold

integer function user_exceeds_threshold(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag, is_ghost, refine

    refine = 0

    if (qval .gt. threshold .and. qval .lt. 1-threshold) then
        refine = 1
    endif

    user_exceeds_threshold = refine

end function user_exceeds_threshold
