SUBROUTINE clawpack5_tag4refinement(mx,my,mbc,meqn, &
     xlower,ylower,dx,dy,blockno,q, &
     refine_threshold, init_flag, &
     tag_for_refinement)
  IMPLICIT NONE

  INTEGER mx,my, mbc, meqn, tag_for_refinement, init_flag
  INTEGER blockno
  DOUBLE PRECISION xlower, ylower, dx, dy
  DOUBLE PRECISION refine_threshold
  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j, mq
  DOUBLE PRECISION xc,yc, qmin, qmax

  tag_for_refinement = 0

  !!  # Refine based only on first variable in system.
  DO mq = 1,meqn
     qmin = 100.d0
     qmax = -100.d0
     DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
           IF (init_flag .EQ. 1) THEN
              xc = xlower + (i-0.5)*dx
              yc = ylower + (j-0.5)*dy
              IF (ABS(xc-0.5) .LT. dy) THEN
                 tag_for_refinement = 1
                 RETURN
              ENDIF
           ELSE
              !! # EXIT immediately IF the refinement criteria is met
              qmin = MIN(q(mq,i,j),qmin)
              qmax = MAX(q(mq,i,j),qmax)
              IF (qmax - qmin .GT. refine_threshold) THEN
                 tag_for_refinement = 1
                 RETURN
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE clawpack5_tag4refinement
