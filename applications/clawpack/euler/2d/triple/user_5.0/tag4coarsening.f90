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

  INTEGER mq
  DOUBLE PRECISION qmin, qmax, dq

  tag_patch = 0
  dq = 0
  DO mq = 1,1
     CALL clawpack5_get_minmax(mx,my,mbc,meqn,mq,q0,dq)
     CALL clawpack5_get_minmax(mx,my,mbc,meqn,mq,q1,dq)
     CALL clawpack5_get_minmax(mx,my,mbc,meqn,mq,q2,dq)
     CALL clawpack5_get_minmax(mx,my,mbc,meqn,mq,q3,dq)
  ENDDO
  IF (dq .LT. coarsen_threshold) THEN
     tag_patch = 1
     RETURN
  ENDIF

END SUBROUTINE clawpack5_tag4coarsening

SUBROUTINE clawpack5_get_minmax(mx,my,mbc,meqn,mq,q,dq)

  IMPLICIT NONE
  INTEGER mx,my,mbc,meqn,mq
  DOUBLE PRECISION qmin,qmax,dq
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  INTEGER i,j
  DOUBLE PRECISION dqi, dqj

  DO i = 1,mx
     DO j = 1,my
        dqi = dabs(q(mq,i+1,j) - q(mq,i-1,j))
        dqj = dabs(q(mq,i,j+1) - q(mq,i,j-1))
        dq  = dmax1(dq,dqi, dqj)
     ENDDO
  ENDDO

END SUBROUTINE clawpack5_get_minmax
