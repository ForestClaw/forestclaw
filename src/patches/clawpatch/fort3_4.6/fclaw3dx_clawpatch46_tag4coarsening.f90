subroutine fclaw3dx_clawpatch46_fort_tag4coarsening(mx,my,mz, mbc,meqn, & 
           xlower,ylower,zlower,dx,dy, dz, blockno, q0, q1, q2, q3, & 
          coarsen_threshold, initflag, tag_patch)
    implicit none

    integer :: mx,my, mz, mbc, meqn, tag_patch, initflag
    integer :: glob, blockno
    double precision :: xlower(0:3), ylower(0:3), zlower, dx, dy, dz
    double precision :: coarsen_threshold
    double precision :: q0(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q1(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q2(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: q3(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: mq
    double precision :: qmin, qmax

    !! # Don't coarsen when initializing the mesh
    if (initflag .ne. 0) then
        tag_patch = 0
        return
    endif

    !! # Assume that we will coarsen a family unless we find a grid
    !! # that doesn't pass the coarsening test.
    tag_patch = 1
    mq = 1
    qmin = q0(1,1,1,mq)
    qmax = q0(1,1,1,mq)

    call fclaw3dx_clawpatch46_test_refine3(glob,blockno,mx,my,mz,mbc,meqn, & 
           mq,q0,qmin,qmax, dx,dy,dz,xlower(0), ylower(0),zlower, &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3dx_clawpatch46_test_refine3(glob,blockno,mx,my,mz,mbc,meqn, & 
           mq,q1,qmin,qmax,dx,dy,dz,xlower(1), ylower(1),  zlower, & 
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3dx_clawpatch46_test_refine3(glob,blockno,mx,my,mz,mbc,meqn, & 
           mq,q2,qmin,qmax,dx,dy,dz,xlower(2), ylower(2),zlower, &
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fclaw3dx_clawpatch46_test_refine3(glob,blockno,mx,my,mz,mbc,meqn, &
           mq,q3,qmin,qmax,dx,dy,dz,xlower(3), ylower(3),zlower, &
           coarsen_threshold,initflag, tag_patch)

end subroutine fclaw3dx_clawpatch46_fort_tag4coarsening

subroutine fclaw3dx_clawpatch46_test_refine3(glob,blockno,mx,my,mz,mbc, & 
      meqn,mq,q, qmin,qmax,dx,dy,dz,xlower,ylower,zlower, &
      coarsen_threshold,init_flag,tag_patch)

    implicit none
    integer :: mx,my,mz,mbc,meqn,mq,tag_patch, init_flag, glob, blockno
    double precision :: coarsen_threshold
    double precision :: qmin,qmax, dx, dy, dz, xlower, ylower,zlower
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    double precision :: xc,yc,zc, quad(-1:1,-1:1,-1:1),qval

    integer i,j, k, ii, jj,kk

    integer :: exceeds_th, fclaw3dx_clawpatch_exceeds_threshold
    logical(kind=4) :: is_ghost, clawpatch3_is_ghost


    do k = 1,mz
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc                
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                zc = zlower + (k-0.5)*dz
                qmin = min(q(i,j,k,mq),qmin)
                qmax = max(q(i,j,k,mq),qmax)
                qval = q(i,j,k,mq)
                is_ghost = clawpatch3_is_ghost(i,j,k,mx,my,mz)
                if (.not. is_ghost) then
                    do ii = -1,1               
                        do jj = -1,1
                            do kk = -1,1
                                quad(ii,jj,kk) = q(i+ii,j+jj,k+kk,mq)
                            end do
                        end do
                    end do
                endif
                exceeds_th = fclaw3dx_clawpatch_exceeds_threshold( & 
                    glob,blockno, qval,qmin,qmax,quad, dx,dy,dz,xc,yc,zc,  &
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

end subroutine fclaw3dx_clawpatch46_test_refine3
