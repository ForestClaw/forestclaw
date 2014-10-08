!!
!!
!!


PROGRAM mix_data
  IMPLICIT NONE

  INTEGER :: lx, ly, iframe

  INTEGER :: i,j, m, mq
  INTEGER :: ngrid, level, mx, my, meqn
  INTEGER :: igrid, ngrids, blocknumber
  DOUBLE PRECISION :: xlow, ylow, dx, dy, t

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: q

  CHARACTER(100) :: fname1, fname2, fname3
  CHARACTER :: c
  INTEGER :: nstp, ipos, idigit, l1, l2

!! Get iframe and file name to create
  OPEN(10,file='diag.dat')
  READ(10,*) iframe
  CLOSE(10)

  fname1 = 'fort.qxxxx'
  fname2 = 'fort.txxxx'
  nstp = iframe
  DO ipos = 10, 7, -1
     idigit = MOD(nstp,10)
     fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     fname2(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     nstp = nstp / 10
  ENDDO

!!   l1 = LEN_TRIM(fname3)
!!   fname3(l1+1:l1+1) = '.'
!!   nstp = iframe
!!   IF (iframe .EQ. 0) THEN
!!      l2 = 1
!!   ELSE
!!      l2 = 0
!!      DO WHILE (nstp .GT. 0)
!!         nstp = nstp/10
!!         l2 = l2 + 1
!!      END DO
!!   ENDIF
!!   nstp = iframe
!!   DO ipos = l1+l2 + 1,l1 + 2, -1
!!      idigit = MOD(nstp,10)
!!      fname3(ipos:ipos) = CHAR(ICHAR('0') + idigit)
!!      nstp = nstp / 10
!!   ENDDO


  OPEN(10,file=fname2)
  READ(10,*) t
  READ(10,*) meqn
  READ(10,*) ngrids
  CLOSE(10)

  OPEN(10,file=fname1)
  OPEN(20,file='mix.out')
  DO igrid = 1, ngrids
     OPEN(10,file=fname1)
     READ(10,*) ngrid
     READ(10,*) level
     READ(10,*) blocknumber
     READ(10,*) mx
     READ(10,*) my
     READ(10,*) xlow
     READ(10,*) ylow
     READ(10,*) dx
     READ(10,*) dy

     ALLOCATE(q(0:mx+1,0:my+1,meqn))


     DO j = 1,my
        DO i = 1,mx
           READ(10,*) (q(i,j,m),m = 1,meqn)
           !! Store third entry as area element, not as a capacity
           q(i,j,3) = dx*dy*q(i,j,3)
        ENDDO
     ENDDO

     !! fname1(17:17) = CHAR(ICHAR('0') + iframe)
     !! fname2(17:17) = CHAR(ICHAR('0') + iframe)

     DO i = 1,mx
        DO j = 1,my
           WRITE(20,120) (q(i,j,m),m=1,meqn)
        ENDDO
        WRITE(20,*) ' '
     ENDDO
     DEALLOCATE(q)
  END DO
  CLOSE(20)
  close(10)
120  FORMAT(3F24.16)

100 FORMAT(A,E10.4)

END PROGRAM mix_data
