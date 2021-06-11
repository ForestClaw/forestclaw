SUBROUTINE clawpack5_tag4coarsening(mx,my,mbc,meqn, &
     xlower,ylower,dx,dy, blockno, q0, q1, q2, q3, &
     coarsen_threshold, init_flag, tag_patch)
  IMPLICIT NONE

  Integer Mx,My, Mbc, Meqn, Tag_Patch
  INTEGER blockno, init_flag
  DOUBLE PRECISION xlower(0:3), ylower(0:3), dx, dy
  DOUBLE PRECISION coarsen_threshold
  DOUBLE PRECISION q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j, mq
  DOUBLE PRECISION qmin, qmax

  !! # Assume that we will coarsen a family unless we find a grid
  !! # that doesn't pass the coarsening test.
  if (init_flag .ne. 0) then
      tag_patch=0
      return
  endif


  tag_patch = 1
  mq = 1
  qmin = q0(mq,1,1)
  qmax = q0(mq,1,1)

  !! # If we find that (qmax-qmin > coarsen_threshold) on any
  !! # grid, we return immediately, since the family will then
  !! # not be coarsened.

  CALL user5_tag_sibling(mx,my,mbc,meqn,mq,q0,qmin,qmax, &
       coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(mx,my,mbc,meqn,mq,q1,qmin,qmax, &
       coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(mx,my,mbc,meqn,mq,q2,qmin,qmax, &
       coarsen_threshold,tag_patch)
  IF (tag_patch == 0) RETURN

  CALL user5_tag_sibling(mx,my,mbc,meqn,mq,q3,qmin,qmax, &
       coarsen_threshold,tag_patch)

END SUBROUTINE clawpack5_tag4coarsening

SUBROUTINE user5_tag_sibling(mx,my,mbc,meqn,mq,q, &
     qmin,qmax,coarsen_threshold,tag_patch)

  IMPLICIT NONE
  INTEGER mx,my,mbc,meqn,mq,tag_patch
  DOUBLE PRECISION coarsen_threshold
  DOUBLE PRECISION qmin,qmax
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  INTEGER i,j

  DO i = 1,mx
     DO j = 1,my
        qmin = min(q(mq,i,j),qmin)
        qmax = max(q(mq,i,j),qmax)
        IF (qmax - qmin.GT. coarsen_threshold) THEN
           !! # We won't coarsen this family because at least one
           !! # grid fails the coarsening test.
           tag_patch = 0
           RETURN
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE user5_tag_sibling
