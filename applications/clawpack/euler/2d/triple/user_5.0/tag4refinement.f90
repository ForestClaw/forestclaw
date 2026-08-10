SUBROUTINE clawpack5_tag4refinement(mx,my,mbc, &
     meqn, xlower,ylower,dx,dy,blockno, &
     q, tag_threshold, init_flag,tag_patch)
  IMPLICIT NONE

  INTEGER mx,my, mbc, meqn, tag_patch, init_flag
  INTEGER blockno
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION tag_threshold
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j, mq
  DOUBLE PRECISION qmin, qmax
  DOUBLE PRECISION dq, dqi, dqj

  tag_patch = 0

  !! # Refine based only on first variable in system.
  qmin = q(1,1,1)
  qmax = q(1,1,1)
  DO j = 2-mbc,my+mbc-1
     DO i = 2-mbc,mx+mbc-1
        dq = 0
        DO mq = 1,1
           dqi = dabs(q(mq,i+1,j) - q(mq,i-1,j))
           dqj = dabs(q(mq,i,j+1) - q(mq,i,j-1))
           dq  = dmax1(dq, dqi, dqj)
           IF (dq .GT. tag_threshold) THEN
              tag_patch = 1
              RETURN
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE clawpack5_tag4refinement
