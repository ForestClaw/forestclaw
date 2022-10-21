!! # check to see if value exceeds threshold

integer function fclaw3dx_user_exceeds_th(blockno,& 
        qval,qmin,qmax,quad, & 
        dx,dy,dz,xc,yc,zc,threshold, &
        init_flag, is_ghost)
    use setprob_mod, only : pi, mapping, init_choice, x0, y0, z0, r0, maxelev
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1,-1:1)
    double precision :: dx,dy, dz,xc, yc,zc
    integer :: blockno, init_flag, refine
    logical(kind=4) :: is_ghost

    refine = 0

    !! Unfortuntely, we are not yet passing in all q values (only q(mq)). 
    !! This is easy to fix. 
    if (.not. is_ghost) then
        if (mapping .eq. 0) then
            if (xc .lt. -0.5) then
                refine = 1
            endif
        elseif (mapping .eq. 1) then
            if (xc .lt. 0.5) then
                refine = 1
            endif
        endif
    endif

    fclaw3dx_user_exceeds_th = refine

end function fclaw3dx_user_exceeds_th
