SUBROUTINE torus5_fort_write_file(matname1, &
     mx,my,meqn,mbc, xlower,ylower, dx,dy, &
     q,error,t,patch_num,level,blockno,mpirank)

  IMPLICIT NONE

  CHARACTER*10 matname1
  INTEGER meqn,mbc,mx,my
  INTEGER patch_num, level, blockno, mpirank
  DOUBLE PRECISION xlower, ylower,dx,dy,t
  DOUBLE PRECISION xc,yc,qc,qexact

  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER matunit1
  INTEGER i,j,mq

  matunit1 = 10
  OPEN(matunit1,file=matname1,position='append');

!!  CALL fc2d_clawpack5_fort_write_grid_header(matunit1, &
!!       mx,my,xlower,ylower, dx,dy,patch_num,level, &
!!       blockno,mpirank)

  IF (meqn .GT. 5) THEN
     WRITE(6,'(A,A,A)') &
          'Warning (fclaw2d_fort_write_grid_header.f) ', &
          ': meqn > 5; change format statement 120.'
     STOP
  ENDIF

  DO j = 1,my
     DO i = 1,mx
        DO mq = 1,meqn
           IF (ABS(q(mq,i,j)) .LT. 1d-99) THEN
              q(mq,i,j) = 0.d0
           ENDIF
        ENDDO
        xc = xlower + (i-0.5)*dx
        yc = ylower + (j-0.5)*dy
        qc = qexact(blockno,xc, yc,t)
        IF (ABS(qc) .LT. 1d-99) THEN
           qc = 0.d0
        ENDIF
        IF (ABS(error(1,i,j)) .LT. 1d-99) THEN
           error(1,i,j) = 0.d0
        ENDIF
        WRITE(matunit1,120) (q(mq,i,j),mq=1,meqn), &
             error(1,i,j), qc
     ENDDO
     WRITE(matunit1,*) ' '
  ENDDO
  !! # This statement is checked above (meqn <= 5)
120 FORMAT (5E26.16)

  CLOSE(matunit1)

END SUBROUTINE torus5_fort_write_file
