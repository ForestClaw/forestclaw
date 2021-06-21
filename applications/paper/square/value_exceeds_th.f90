!! # check to see if value exceeds threshold

logical function value_exceeds_th(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshhold)
    implicit none
    
    double precision :: qval,qmin,qmax,threshhold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno

    integer initchoice
    common /initchoice_comm/ initchoice

    integer example
    common /example_comm/ example      

    integer :: mq
    logical :: refine

    refine = .false.

    mq = 1
    if (initchoice .le. 1) then
        if (example .eq. 0) then
            refine = qval .gt. threshhold .and. &
                     qval .lt. 1-threshhold
        else
            refine = qval .gt. threshhold
        endif
    elseif (initchoice .eq. 2) then
        refine = qval .gt.  threshhold              
    else                
        write(6,'(A,A)') 'Refining not yet defined for example > 0'
        stop
    endif


    value_exceeds_th = refine

end function value_exceeds_th

