!!
!! gfortran -o compute_filament -fdefault-real-8 -fdefault-double-8 diag_new.f90 compute_filament.f90
!!


PROGRAM compute_filament
    USE diag
    IMPLICIT NONE

    INTEGER :: iframe

    INTEGER :: i,j, k,kmax
    INTEGER :: ngrid, level, mx, my, meqn, ngrids
    integer :: blockno, mpirank
    DOUBLE PRECISION :: xlow, ylow, dx, dy, t

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f0, f1, f2, dA
    DOUBLE PRECISION s, err

    CHARACTER(100) :: fname1, fname2
    CHARACTER :: c
    INTEGER :: nstp, ipos, idigit, l1, l2
    DOUBLE PRECISION :: diagout(3), fila_t0(100)
    LOGICAL :: linit

    INTEGER m, igrid
    double precision t0, t1

    iframe = 0
    fname1 = 'fort.qxxxx'
    fname2 = 'fort.txxxx'
    nstp = iframe
    DO ipos = 10, 7, -1
       idigit = MOD(nstp,10)
       fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
       fname2(ipos:ipos) = CHAR(ICHAR('0') + idigit)
       nstp = nstp / 10
    END DO

    OPEN(10,file=fname2)
    READ(10,*) t0
    READ(10,*) meqn
    read(10,*) ngrids
    CLOSE(10)

!! Get mx,my
    OPEN(10,file=fname1)  
    READ(10,*) ngrid
    READ(10,*) level
    READ(10,*) blockno
    READ(10,*) mpirank
    READ(10,*) mx
    READ(10,*) my
    CLOSE(10)


    kmax = ngrids*mx*my
    ALLOCATE(f0(kmax), f1(kmax),f2(kmax),dA(kmax))

    k = 1
    s = 0
    OPEN(10,file=fname1)
    DO igrid = 1, ngrids
        READ(10,*) ngrid
        READ(10,*) level
        READ(10,*) blockno
        READ(10,*) mpirank
        READ(10,*) mx
        READ(10,*) my
        READ(10,*) xlow
        READ(10,*) ylow
        READ(10,*) dx
        READ(10,*) dy

        DO j = 1,my
            DO i = 1,mx
                READ(10,*) f0(k),f2(k),err, dA(k)
                s = s + f0(k)*dA(k)
                k = k + 1
            END DO
        END DO

    END DO
100 FORMAT(A,E16.8)
    CLOSE(10)

    write(6,101) 'Mass at time ',t0,' : ', s
101 format(A,F6.4,A,F24.16)    
    linit = .TRUE.
    !! This call sets fila_t0
    CALL filament_diag(Kmax,f0,dA,fila_t0,linit)

!! -----------------------------------------------------------
!! Now read data at T/2
!! -----------------------------------------------------------
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
    END DO

    OPEN(10,file=fname2)
    READ(10,*) t1
    READ(10,*) meqn
    read(10,*) ngrids
    CLOSE(10)

    !! For adaptive runs, the number of grids at time 'iframe > ' may be different
    !! than the number of grids at time iframe=0
    DEALLOCATE(f0, f1,f2,dA)

    kmax = ngrids*mx*my
    ALLOCATE(f0(kmax), f1(kmax),f2(kmax),dA(kmax))

    s = 0
    k = 1
    OPEN(10,file=fname1)
    DO igrid = 1,ngrids
        READ(10,*) ngrid
        READ(10,*) level
        read(10,*) blockno
        read(10,*) mpirank
        READ(10,*) mx
        READ(10,*) my
        READ(10,*) xlow
        READ(10,*) ylow
        READ(10,*) dx
        READ(10,*) dy

        DO j = 1,my
            DO i = 1,mx
                !! Read data from fort.q file
                !! f1(k) : Computed solution
                !! f2(k) : Exact solution (or computed solution, if error is not computed)
                !! err   : Error 
                !! dA(k) : Area element
                READ(10,*) f1(k),f2(k),err, dA(k)
                s = s + f1(k)*dA(k)
                k = k + 1
            END DO
        END DO
    END DO
    CLOSE(10)

    write(6,101) 'Mass at time ',t1,' : ', s
    linit = .FALSE.
    !! This call uses initial data in fila_t0
    CALL filament_diag(Kmax,f1,dA,fila_t0,linit)

END PROGRAM compute_filament
