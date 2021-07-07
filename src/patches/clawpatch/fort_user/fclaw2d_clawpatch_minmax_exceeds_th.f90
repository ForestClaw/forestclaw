!! # check to see if value exceeds threshold

logical function fclaw2d_clawpatch_minmax_exceeds_th(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, is_ghost, init_flag

    logical :: refine

    refine = .false.

    if (qmax-qmin .gt. threshold) then
        refine = .true.
    endif

    fclaw2d_clawpatch_minmax_exceeds_th = refine

end function fclaw2d_clawpatch_minmax_exceeds_th
