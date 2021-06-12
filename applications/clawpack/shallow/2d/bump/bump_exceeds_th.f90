!! # check to see if value exceeds threshold

logical function bump_exceeds_th(qval,qmin,qmax,threshhold)
    implicit none
    
    double precision qval,qmin,qmax,threshhold
    logical refine

    refine = .false.
    if (qval .gt. threshhold) then
        refine = .true.
    endif

    bump_exceeds_th = refine

end function bump_exceeds_th