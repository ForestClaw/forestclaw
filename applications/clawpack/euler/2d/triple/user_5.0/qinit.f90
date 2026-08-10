SUBROUTINE clawpack5_qinit(meqn,mbc, &
     mx,my,xlower,ylower,dx,dy,q,maux,aux)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy
  DOUBLE PRECISION   q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  INTEGER blockno, fc2d_clawpack5_get_block

  DOUBLE PRECISION qin(5),qout(5),rinf,vinf,einf
  COMMON /comic/ qin,qout
  COMMON /cominf/ rinf,vinf,einf

  INTEGER i,j,m
  DOUBLE PRECISION xclow,yclow,win

  blockno = fc2d_clawpack5_get_block()

  DO i = 1-mbc,mx+mbc
     xclow = xlower + (i-1)*dx
     DO j = 1-mbc,my+mbc
        yclow = ylower + (j-1)*dy
        CALL cellave2(blockno,xclow,yclow,dx,dy,win)
        DO m = 1,meqn
           q(m,i,j) = win*qin(m) + (1.d0-win)*qout(m)
        ENDDO
     ENDDO

     IF (xclow .LT. 0.2d0) THEN
        DO j=1-mbc,my+mbc
           q(1,i,j) = rinf
           q(2,i,j) = rinf*vinf
           q(3,i,j) = 0.d0
           q(4,i,j) = einf
           IF (meqn .EQ. 5) THEN
              q(5,i,j) = 1.d0
           ENDIF
        ENDDO
     END IF

     IF (meqn .EQ. 5) THEN
        IF (xclow .LT. 0.5d0) THEN
           DO j = 1,my
              q(5,i,j) = 2.d0*q(5,i,j)
           ENDDO
        END IF
     END IF
  ENDDO

  RETURN
END SUBROUTINE clawpack5_qinit
