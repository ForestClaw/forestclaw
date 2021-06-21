SUBROUTINE clawpatch5_tag4coarsening(mx,my,mbc,meqn, &
     xlower,ylower,dx,dy, blockno, q0, q1, q2, q3, &
     coarsen_threshold, init_flag, tag_patch)
  IMPLICIT NONE

  INTEGER :: mx,my, mbc, meqn, tag_patch
  INTEGER :: blockno, init_flag
  DOUBLE PRECISION :: xlower(0:3), ylower(0:3), dx, dy
  DOUBLE PRECISION :: coarsen_threshold
  DOUBLE PRECISION :: q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION :: q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION :: q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION :: q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER :: mq
  DOUBLE PRECISION :: qmin, qmax

  if (init_flag .ne. 0) then
      tag_patch = 0
      return
  endif

  !! # Assume that we will coarsen a family unless we find a grid
  !! # that doesn't pass the coarsening test.
  tag_patch = 1
  mq = 1
  qmin = q0(mq,1,1)
  qmax = q0(mq,1,1)

  !! # If we find that (qmax-qmin > coarsen_threshold) on any
  !! # grid, we return immediately, since the family will then
  !! # not be coarsened.

  CALL user5_tag_sibling(blockno, mx,my,mbc,meqn,mq,q0,qmin,qmax, &
       xlower(0), ylower(0), dx, dy, coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(blockno, mx,my,mbc,meqn,mq,q1,qmin,qmax, &
       xlower(1), ylower(1), dx, dy, coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(blockno, mx,my,mbc,meqn,mq,q2,qmin,qmax, &
       xlower(2), ylower(2), dx, dy, coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(blockno, mx,my,mbc,meqn,mq,q3,qmin,qmax, &
       xlower(3), ylower(3), dx, dy, coarsen_threshold,tag_patch)

END SUBROUTINE clawpatch5_tag4coarsening

SUBROUTINE user5_tag_sibling(blockno, mx,my,mbc,meqn,mq,q, &
     qmin,qmax,xlower, ylower, dx, dy, & 
     coarsen_threshold,tag_patch)

    IMPLICIT NONE
    INTEGER :: mx,my,mbc,meqn,mq,tag_patch,blockno
    double precision :: dx,dy,xlower, ylower
    DOUBLE PRECISION :: coarsen_threshold
    DOUBLE PRECISION :: qmin,qmax
    DOUBLE PRECISION :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    DOUBLE PRECISION :: quad(-1:1,-1:1), xc, yc
    LOGICAL :: exceeds_th, gradient_exceeds_th
    INTEGER :: i,j,ii,jj

    DO i = 1,mx
        DO j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(mq, i,j),qmin)
            qmax = max(q(mq, i,j),qmax)
            do ii = -1,1
               do jj = -1,1
                  quad(ii,jj) = q(mq,i+ii,j+jj)
               end do
            end do

            exceeds_th = gradient_exceeds_th(blockno,  & 
                  q(mq, i,j),qmin,qmax,quad,dx,dy,xc,yc,  & 
                  coarsen_threshold)
            if (exceeds_th) then
!!              # We won't coarsen this family because at least one
!!              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif

        END DO
    END DO

END SUBROUTINE user5_tag_sibling
