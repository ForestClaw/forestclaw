SUBROUTINE clawpatch5_tag4refinement(mx,my,mbc,meqn, &
           xlower,ylower,dx,dy,blockno, q,refine_threshold, &
           init_flag, tag_patch)
    IMPLICIT NONE

    INTEGER :: mx,my, mbc, meqn, tag_patch, init_flag
    INTEGER :: blockno
    DOUBLE PRECISION :: xlower, ylower, dx, dy
    DOUBLE PRECISION :: refine_threshold
    DOUBLE PRECISION :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER :: i,j, mq,m, ii, jj
    DOUBLE PRECISION :: xc,yc, qmin, qmax, quad(-1:1,-1:1)
    logical :: exceeds_th, gradient_exceeds_th

    tag_patch = 0

    !! # Refined based on dq/dr (radial derivative)

    qmin = q(mq,1,1)
    qmax = q(mq,1,1)
    mq = 1
    DO i = 1,mx
       DO j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            do ii = -1,1
                do jj = -1,1
                    quad(ii,jj) = q(mq,i+ii,j+jj)
                end do
            end do
            exceeds_th = gradient_exceeds_th(blockno,  &                                       
                    q(mq,i,j),qmin,qmax,quad, dx,dy,xc,yc, &
                    refine_threshold)
            if (exceeds_th) then
                tag_patch = 1
                return
            endif
        end do
    END DO

END SUBROUTINE clawpatch5_tag4refinement
