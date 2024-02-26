subroutine fclaw3dx_clawpatch46_fort_tag4coarsening(mx,my,mz, mbc,meqn, & 
           xlower,ylower,zlower,dx,dy, dz, blockno, q0, q1, q2, q3, &
           coarsen_threshold, initflag, tag_patch)
    implicit none

    integer :: mx,my, mz, mbc, meqn, tag_patch, initflag
    integer :: blockno
    double precision :: xlower(0:3), ylower(0:3), zlower, dx, dy, dz
    double precision :: coarsen_threshold
    double precision :: q0(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q1(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q2(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q3(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)


    integer :: mq
    double precision :: qmin(meqn), qmax(meqn)

    !! # Don't coarsen when initializing the mesh
    if (initflag .ne. 0) then
        tag_patch = 0
        return
    endif

    !! # Assume that we will coarsen a family unless we find a grid
    !! # that doesn't pass the coarsening test.
    tag_patch = 1
    do mq = 1,meqn
        qmin(mq) = q0(1,1,1,mq)
        qmax(mq) = q0(1,1,1,mq)
    end do

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q0,qmin,qmax, dx,dy,dz,xlower(0), ylower(0),zlower,  &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q1,qmin,qmax,dx,dy,dz,xlower(1), ylower(1),  zlower, & 
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q2,qmin,qmax,dx,dy,dz,xlower(2), ylower(2),zlower,  &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q3,qmin,qmax,dx,dy,dz,xlower(3), ylower(3),zlower, &
           coarsen_threshold,initflag, tag_patch)

end subroutine fclaw3dx_clawpatch46_fort_tag4coarsening

subroutine fclaw3d_clawpatch46_fort_tag4coarsening(mx,my,mz, mbc,meqn, & 
           xlower,ylower,zlower,dx,dy, dz, blockno, &
           q0, q1, q2, q3, q4, q5, q6, q7, &
           coarsen_threshold, initflag, tag_patch)
    implicit none

    integer :: mx,my, mz, mbc, meqn, tag_patch, initflag
    integer :: blockno
    double precision :: xlower(0:7), ylower(0:7), zlower(0:7), dx, dy, dz
    double precision :: coarsen_threshold
    double precision :: q0(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q1(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q2(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q3(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q4(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q5(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q6(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q7(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)


    integer :: mq
    double precision :: qmin(meqn), qmax(meqn)

    !! # Don't coarsen when initializing the mesh
    if (initflag .ne. 0) then
        tag_patch = 0
        return
    endif

    !! # Assume that we will coarsen a family unless we find a grid
    !! # that doesn't pass the coarsening test.
    tag_patch = 1
    do mq = 1,meqn
        qmin(mq) = q0(1,1,1,mq)
        qmax(mq) = q0(1,1,1,mq)
    end do

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q0,qmin,qmax, dx,dy,dz,xlower(0), ylower(0),zlower(0),  &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q1,qmin,qmax,dx,dy,dz,xlower(1), ylower(1), zlower(1), & 
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, & 
           q2,qmin,qmax,dx,dy,dz,xlower(2), ylower(2),zlower(2),  &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q3,qmin,qmax,dx,dy,dz,xlower(3), ylower(3),zlower(3), &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q4,qmin,qmax,dx,dy,dz,xlower(4), ylower(4),zlower(4), &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q5,qmin,qmax,dx,dy,dz,xlower(5), ylower(5),zlower(5), &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q6,qmin,qmax,dx,dy,dz,xlower(6), ylower(6),zlower(6), &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc,meqn, &
           q7,qmin,qmax,dx,dy,dz,xlower(7), ylower(7),zlower(7), &
           coarsen_threshold,initflag, tag_patch)

end subroutine fclaw3d_clawpatch46_fort_tag4coarsening

subroutine fclaw3d_clawpatch46_test_refine3(blockno,mx,my,mz,mbc, & 
      meqn,q, qmin,qmax,dx,dy,dz,xlower,ylower,zlower, &
      coarsen_threshold,init_flag,tag_patch)

    implicit none
    integer :: mx,my,mz,mbc,meqn,tag_patch, init_flag, blockno
    double precision :: coarsen_threshold
    double precision :: dx, dy, dz, xlower, ylower,zlower
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: qmin(meqn),qmax(meqn)

    !! Dummy variables
    double precision :: xc,yc,zc, quad(-1:1,-1:1,-1:1,meqn),qval(meqn)

    integer i,j, k, ii, jj,kk, mq

    integer :: exceeds_th, fclaw3d_clawpatch_tag_criteria
    logical(kind=4) :: is_ghost, clawpatch3_is_ghost


    do k = 1-mbc,mz+mbc
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc                
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                zc = zlower + (k-0.5)*dz
                do mq = 1,meqn
                    qval(mq) = q(i,j,k,mq)
                    qmin(mq) = min(q(i,j,k,mq),qmin(mq))
                    qmax(mq) = max(q(i,j,k,mq),qmax(mq))
                end do
                is_ghost = clawpatch3_is_ghost(i,j,k,mx,my,mz)
                if (.not. is_ghost) then
                    do ii = -1,1               
                        do jj = -1,1
                            do kk = -1,1
                                do mq = 1,meqn
                                    quad(ii,jj,kk,mq) = q(i+ii,j+jj,k+kk,mq)
                                end do
                            end do
                        end do
                    end do
                endif
                exceeds_th = fclaw3d_clawpatch_tag_criteria( & 
                    blockno, qval,qmin,qmax,quad, dx,dy,dz,xc,yc,zc,  &
                    coarsen_threshold, init_flag, is_ghost)
            
                !! # -1 : Not conclusive (possibly ghost cell) (do not tag for coarsening)
                !! # 0  : Does not pass threshold (tag for coarsening)      
                !! # 1  : Passes threshold (do not tag for coarsening)
                !! # Note : exceeds_th = -1 leads to over-refining, so it is 
                !! # ignored here.  Logic of regridding (coarsening then 
                !! # refining) isn't clear.
                if (exceeds_th .gt. 0) then
                    tag_patch = 0
                    return
                endif
            end do
        end do
    end do

end subroutine fclaw3d_clawpatch46_test_refine3
