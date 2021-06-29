!! # check to see if value exceeds threshold

logical function minmax_exceeds_th(blockno,& 
                                   qval,qmin,qmax,quad, & 
                                   dx,dy,xc,yc,threshhold)
    implicit none
    
    double precision :: qval,qmin,qmax,threshhold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno

    logical :: refine

    refine = .false.
    if (qmax-qmin .gt. threshhold) then
        refine = .true.
    endif

    minmax_exceeds_th = refine

end function minmax_exceeds_th