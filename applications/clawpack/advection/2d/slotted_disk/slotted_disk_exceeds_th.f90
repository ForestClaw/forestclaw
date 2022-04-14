!! # check to see if value exceeds threshold

integer function user_exceeds_threshold(glob, blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: glob, blockno, init_flag, refine
    logical(kind=4) :: is_ghost

    refine = 0

    if (.not. is_ghost) then
        if (qval .gt. threshold .and. qval .lt. 1-threshold) then
            refine = 1
        endif
    endif

    user_exceeds_threshold = refine

end function user_exceeds_threshold
