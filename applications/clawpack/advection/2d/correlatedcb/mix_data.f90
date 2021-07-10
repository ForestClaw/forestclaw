!!
!! gfortran -o mix_data -fdefault-real-8 -fdefault-double-8 -g mix_data.f90
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

  INTEGER :: mpirank

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
     READ(10,*) mpirank
     READ(10,*) mx
     READ(10,*) my
     READ(10,*) xlow
     READ(10,*) ylow
     READ(10,*) dx
     READ(10,*) dy

     ALLOCATE(q(mx,my,meqn))

     DO j = 1,my
        DO i = 1,mx
           READ(10,*) (q(i,j,m),m = 1,meqn)
        ENDDO
     ENDDO

     DO i = 1,mx
        DO j = 1,my
           WRITE(20,120) (q(i,j,m),m=1,meqn)
        ENDDO
        WRITE(20,*) ' '
     ENDDO
     DEALLOCATE(q)
  END DO
  CLOSE(20)
  CLOSE(10)
120  FORMAT(3E24.16)

100 FORMAT(A,E10.4)

END PROGRAM mix_data
