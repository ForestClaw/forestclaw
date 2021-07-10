!!
!! gfortran -o compute_mixing -fdefault-real-8 -fdefault-double-8 diag_new.f90 compute_mixing.f90
!!

PROGRAM compute_diag
  USE diag
  IMPLICIT NONE

  INTEGER :: iframe

  INTEGER :: i,j, k,kmax, block_number
  INTEGER :: ngrid, level, mx, my, meqn, ngrids, igrid
  DOUBLE PRECISION :: xlow, ylow, dx, dy, t
  INTEGER :: mpirank

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f1, f2, dA
  DOUBLE PRECISION s

  CHARACTER(100) :: fname1, fname2
  CHARACTER :: c
  INTEGER :: nstp, ipos, idigit, l1, l2
  DOUBLE PRECISION :: diagout(3), fila_t0(100)
  LOGICAL :: linit

  INTEGER m

!! Size of Lat-long grid
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

!! Get mx,my
  OPEN(10,file=fname1)
  READ(10,*) ngrid
  READ(10,*) level
  READ(10,*) block_number
  READ(10,*) mpirank
  READ(10,*) mx
  READ(10,*) my
  CLOSE(10)

  kmax = ngrids*mx*my

  ALLOCATE(f1(kmax),f2(kmax),dA(kmax))

  k = 1
  s = 0
  OPEN(10,file=fname1)
  DO igrid = 1, ngrids
     READ(10,*) ngrid
     READ(10,*) level
     READ(10,*) block_number
     READ(10,*) mpirank
     READ(10,*) mx
     READ(10,*) my
     READ(10,*) xlow
     READ(10,*) ylow
     READ(10,*) dx
     READ(10,*) dy

     DO j = 1,my
        DO i = 1,mx
           READ(10,*) f1(k),f2(k),dA(k)
           s = s + f1(k)*dA(k)
           k = k + 1
!!           IF (k .GT. kmax) THEN
!!              WRITE(6,*) 'k > kmax; kmax = ', kmax
!!              STOP
!!           ENDIF
        ENDDO
     ENDDO
  END DO

100 FORMAT(A,E16.8)
  CLOSE(10)

  CALL correlation_diag(f1,f2,Kmax,dA,diagout)

  OPEN(10,file='diagfile.in')
  WRITE(10,200) my, (diagout(m),m=1,3)
  CLOSE(10)

200 FORMAT(I5,3E16.8)

!!   linit = .TRUE.
!!   CALL filament_diag(Kmax,f1,dA,fila_t0,linit)
!!
!!   linit = .FALSE.
!!   CALL filament_diag(Kmax,f1,dA,fila_t0,linit)

END PROGRAM compute_diag
