!! # check to see if value exceeds threshold

integer function user_exceeds_threshold(blockno,meqn,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold, ivar_threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    integer :: meqn,ivar_threshold
    double precision :: qval(meqn),qmin(meqn),qmax(meqn),threshold
    double precision :: quad(-1:1,-1:1,meqn)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag, refine
    logical(kind=4) :: is_ghost

    integer :: mq

    refine = 0

    mq = ivar_threshold

    if (.not. is_ghost) then
        if (qval(mq) .gt. threshold .and. qval(mq) .lt. 1-threshold) then
            refine = 1
        endif
    endif

    user_exceeds_threshold = refine

end function user_exceeds_threshold
