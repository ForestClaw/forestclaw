!! # check to see if value exceeds threshold

logical function fclaw2d_clawpatch_difference_exceeds_th(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshold,&
                                     init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag
    logical(kind=4) :: is_ghost

    double precision :: dqx, dqy, dq
    logical :: refine

    if (is_ghost) then
!!      # quad likely has uninitialized values. Don't refine
        fclaw2d_clawpatch_difference_exceeds_th = .false.
        return
    endif

    dqx = abs(quad(1,0) - quad(-1,0))
    dqy = abs(quad(0,1) - quad(0,-1))
    dq  = max(dqx, dqy)

    refine = .false.
    if (dq .gt. threshold) then
        refine = .true.
    endif

    fclaw2d_clawpatch_difference_exceeds_th = refine

end function fclaw2d_clawpatch_difference_exceeds_th
