!! # check to see if value exceeds threshold

logical function transport_exceeds_threshold(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,tag_threshold)
    implicit none
    
    double precision :: qval,qmin,qmax,tag_threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno
    logical :: refine

    integer :: refine_criteria
    common /com_refine/ refine_criteria

    logical :: exceeds_th
    logical :: fclaw2d_clawpatch_gradient_exceeds_th
    logical :: fclaw2d_clawpatch_difference_exceeds_th
    logical :: fclaw2d_clawpatch_value_exceeds_th
    logical :: fclaw2d_clawpatch_minmax_exceeds_th

    if (refine_criteria .eq. 0) then 
        exceeds_th = fclaw2d_clawpatch_value_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 1) then
        exceeds_th = fclaw2d_clawpatch_minmax_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 2) then
        exceeds_th = fclaw2d_clawpatch_difference_exceeds_th(blockno, &
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 3) then
        exceeds_th = fclaw2d_clawpatch_gradient_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    end if

    transport_exceeds_threshold = exceeds_th

end function transport_exceeds_threshold
