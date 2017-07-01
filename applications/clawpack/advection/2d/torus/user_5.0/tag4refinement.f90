SUBROUTINE clawpack5_tag4refinement(mx,my,mbc,meqn, &
     xlower,ylower,dx,dy,blockno,q,refine_threshold, &
     init_flag,tag_for_refinement)
  IMPLICIT NONE

  INTEGER mx,my, mbc, meqn, tag_for_refinement, init_flag
  INTEGER blockno
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION refine_threshold
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j, mq
  DOUBLE PRECISION xc,yc, qmin, qmax,qx,qy

  tag_for_refinement = 0

  !! # Refine based only on first variable in system.
  qmin = 0.025d0
  qmax = 0.975d0
  DO mq = 1,meqn
     DO i = 1,mx
        DO j = 1,my
           qx = (q(1,i+1,j)-q(1,i-1,j))/(2*dx)
           qy = (q(1,i,j+1)-q(1,i,j-1))/(2*dy)
           IF (ABS(qx) .GT. refine_threshold .OR. &
                ABS(qy) .GT. refine_threshold) THEN
              tag_for_refinement = 1
              RETURN
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE clawpack5_tag4refinement
