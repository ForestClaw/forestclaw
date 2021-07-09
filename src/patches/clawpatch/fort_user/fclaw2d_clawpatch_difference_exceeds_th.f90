!! # check to see if value exceeds threshold

integer function fclaw2d_clawpatch_difference_exceeds_th(blockno,& 
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
    integer :: refine

    if (is_ghost) then
!!      # quad may have uninitialized values;  test inconclusive
        fclaw2d_clawpatch_difference_exceeds_th = -1
        return
    endif

    dqx = abs(quad(1,0) - quad(-1,0))
    dqy = abs(quad(0,1) - quad(0,-1))
    dq  = max(dqx, dqy)

    refine = 0
    if (dq .gt. threshold) then
        refine = 1
    endif

    fclaw2d_clawpatch_difference_exceeds_th = refine

end function fclaw2d_clawpatch_difference_exceeds_th
