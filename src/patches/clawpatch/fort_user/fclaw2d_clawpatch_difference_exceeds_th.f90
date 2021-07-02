!! # check to see if value exceeds threshold

logical function fclaw2d_clawpatch_difference_exceeds_th(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshhold)
    implicit none
    
    double precision :: qval,qmin,qmax,threshhold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno

    double precision :: dqx, dqy, dq
    logical :: refine

    refine = .false.

    dqx = abs(quad(1,0) - quad(-1,0))
    dqy = abs(quad(0,1) - quad(0,-1))
    dq  = max(dqx, dqy)

    if (dq .gt. threshhold) then
        refine = .true.
    endif

    fclaw2d_clawpatch_difference_exceeds_th = refine

end function fclaw2d_clawpatch_difference_exceeds_th
