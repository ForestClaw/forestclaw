!! # We tag for coarsening if this coarsened patch isn't tagged for refinement
SUBROUTINE clawpack5_tag4coarsening(mx,my,mbc,meqn, &
     xlower,ylower,dx,dy, blockno, q0, q1, q2, q3, &
     coarsen_threshold, tag_patch)
  IMPLICIT NONE

  INTEGER mx,my, mbc, meqn, tag_patch
  INTEGER blockno
  DOUBLE PRECISION xlower(0:3), ylower(0:3), dx, dy
  DOUBLE PRECISION coarsen_threshold
  DOUBLE PRECISION q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j, mq
  DOUBLE PRECISION qmin, qmax

  qmin = 0.025d0
  qmax = 0.975d0
  tag_patch = 1
  CALL torus5_tag_get_minmax(mx,my,mbc,meqn,q0,dx,dy, &
       qmin,qmax,coarsen_threshold,tag_patch)
  IF (tag_patch .EQ. 0) THEN
     RETURN
  ENDIF
  CALL torus5_tag_get_minmax(mx,my,mbc,meqn,q1,dx,dy, &
       qmin,qmax,coarsen_threshold,tag_patch)
  IF (tag_patch .EQ. 0) THEN
     RETURN
  ENDIF
  CALL torus5_tag_get_minmax(mx,my,mbc,meqn,q2,dx,dy, &
       qmin,qmax,coarsen_threshold,tag_patch)
  IF (tag_patch .EQ. 0) THEN
     RETURN
  ENDIF
  CALL torus5_tag_get_minmax(mx,my,mbc,meqn,q3,dx,dy, &
       qmin,qmax,coarsen_threshold,tag_patch)
  IF (tag_patch .EQ. 0) THEN
     RETURN
  ENDIF
  RETURN

END SUBROUTINE clawpack5_tag4coarsening

SUBROUTINE torus5_tag_get_minmax(mx,my,mbc,meqn,q, &
     dx,dy,qmin,qmax,coarsen_threshold,tag_patch)

  IMPLICIT NONE
  INTEGER mx,my,mbc,meqn, tag_patch
  DOUBLE PRECISION dx,dy,qx,qy,coarsen_threshold
  DOUBLE PRECISION qmin,qmax
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  INTEGER i,j,mq

  mq = 1
  DO i = 1,mx
     DO j = 1,my
        qx = (q(1,i+1,j)-q(1,i-1,j))/(2*dx)
        qy = (q(1,i,j+1)-q(1,i,j-1))/(2*dy)
        IF (ABS(qx) .GT. coarsen_threshold .OR. &
             ABS(qy) .GT. coarsen_threshold) THEN
           tag_patch = 0
           RETURN
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE torus5_tag_get_minmax
