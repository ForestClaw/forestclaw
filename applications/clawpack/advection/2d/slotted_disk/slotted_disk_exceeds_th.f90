!! # check to see if value exceeds threshold

logical(kind=4) function user_exceeds_threshold(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag

    logical(kind=4) :: refine, is_ghost

    refine = .false.
    if (qval .gt. threshold .and. qval .lt. 1-threshold) then
        refine = .true.
    endif

    user_exceeds_threshold = refine

end function user_exceeds_threshold
