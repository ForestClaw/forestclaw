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
    double precision :: qmin(meqn), qmax(meqn)

    integer :: exceeds_th, fclaw3dx_clawpatch_tag_criteria
    integer :: ii,jj,kk
    double precision :: xc,yc,zc, quad(-1:1,-1:1,-1:1,meqn), qval(meqn)

    logical(kind=4) :: is_ghost, clawpatch3_is_ghost

    !! # Assume that we won't refine      
    tag_patch = 0

    !! # Default : Refinement based only on first variable in system.  
    !! # Users can modify this by creating a local copy of this routine
    !! # and the corresponding tag4coarsening routine.

    do mq = 1,meqn
        qmin(mq) = q(1,1,1,mq)
        qmax(mq) = q(1,1,1,mq)
    end do
    do k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy 
                zc = zlower + (k-0.5)*dz
                do mq = 1,meqn
                    qval(mq) = q(i,j,k,mq)
                    qmin(mq) = min(qmin(mq),q(i,j,k,mq))                    
                    qmax(mq) = max(qmax(mq),q(i,j,k,mq))
                end do
                is_ghost = clawpatch3_is_ghost(i,j,k, mx,my,mz)
                if (.not. is_ghost) then
                    do jj = -1,1
                        do ii = -1,1
                            do kk = -1,1
                                do mq = 1,meqn
                                    quad(ii,jj,kk,mq) = q(i+ii,j+jj,k+kk,mq)
                                end do
                            end do
                        end do
                    end do
                endif
                exceeds_th = fclaw3dx_clawpatch_tag_criteria(& 
                      blockno, qval,qmin,qmax,quad, dx,dy,dz,xc,yc,zc, & 
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
logical(kind=4) function clawpatch3_is_ghost(i,j,k, mx,my,mz)
    implicit none

    integer :: i, j, k, mx, my, mz
    logical(kind=4) :: interior

    interior = (i .ge. 1 .and. i .le. mx .and. & 
                j .ge. 1 .and. j .le. my .and. &
                k .ge. 1 .and. k .le. mz)
    clawpatch3_is_ghost = .not. interior

    return 

end function clawpatch3_is_ghost


