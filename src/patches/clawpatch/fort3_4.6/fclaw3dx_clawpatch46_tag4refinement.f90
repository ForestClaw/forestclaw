subroutine fclaw3dx_clawpatch46_fort_tag4refinement(mx,my,mz,mbc, & 
    meqn, xlower,ylower,zlower, dx,dy,dz, blockno, & 
    q, tag_threshold, init_flag, tag_patch)
    implicit none

    integer :: mx,my, mz, mbc, meqn, tag_patch, init_flag
    integer :: blockno
    double precision :: xlower, ylower, zlower, dx, dy, dz
    double precision :: tag_threshold
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: i,j,k, mq
    double precision :: qmin, qmax

    integer :: exceeds_th, fclaw3dx_clawpatch_exceeds_threshold
    integer :: ii,jj
    double precision :: xc,yc,quad(-1:1,-1:1), qval

    logical(kind=4) :: is_ghost, clawpatch3_is_ghost

    !! # Assume that we won't refine      
    tag_patch = 0

    !! # Default : Refinement based only on first variable in system.  
    !! # Users can modify this by creating a local copy of this routine
    !! # and the corresponding tag4coarsening routine.

    mq = 1
    qmin = q(1,1,1,mq)
    qmax = q(1,1,1,mq)
    do k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy                
                qmin = min(qmin,q(i,j,k,mq))
                qmax = max(qmax,q(i,j,k,mq))
                qval = q(i,j,k,mq)
                is_ghost = clawpatch3_is_ghost(i,j,mx,my)
                if (.not. is_ghost) then
                    do jj = -1,1
                        do ii = -1,1
                            quad(ii,jj) = q(i+ii,j+jj,k,mq)
                        end do
                    end do
                endif
                exceeds_th = fclaw3dx_clawpatch_exceeds_threshold( & 
                      blockno, qval,qmin,qmax,quad, dx,dy,xc,yc,  & 
                      tag_threshold,init_flag, is_ghost)
                !! # -1 : Not conclusive (possibly ghost cell); don't tag for refinement
                !! # 0  : Does not pass threshold (don't tag for refinement)      
                !! # 1  : Passes threshold (tag for refinement)
                if (exceeds_th .gt. 0) then
                    tag_patch = 1
                    return
                endif
            end do
        end do
    end do
end subroutine  fclaw3dx_clawpatch46_fort_tag4refinement

!! # We may want to check ghost cells for tagging.  
logical(kind=4) function clawpatch3_is_ghost(i,j,mx,my)
    implicit none

    integer :: i, j, mx, my
    logical(kind=4) :: is_ghost

    is_ghost = .false.
    if (i .lt. 1 .or. j .lt. 1) then
        is_ghost = .true.
    elseif (i .gt. mx .or. j .gt. my) then
        is_ghost = .true.
    end if

    clawpatch3_is_ghost = is_ghost

    return 

end function clawpatch3_is_ghost


