!!
!! Compare two fort files in two different directories.
!! To compile :
!!      gfortran -o compare_files compare_files.f90
!!
!! To use : Create a file called

PROGRAM compare_files2d_3dx
    IMPLICIT NONE

    integer :: i,j,k, n, m

    INTEGER :: iframe
    !!  CHARACTER(100) :: dir
    INTEGER :: ngrids, meqn
    INTEGER :: ngrids1, meqn1
    INTEGER :: ngrids2, meqn2
    INTEGER :: ngrid, level, mx, my, mz, blockno, mpirank
    INTEGER :: ngrid1, level1, mx1, my1, mz1, blockno1, mpirank1
    INTEGER :: ngrid2, level2, mx2, my2, mz2, blockno2, mpirank2
    integer :: zlevel

    DOUBLE PRECISION :: xlow, ylow, zlow, dx, dy, dz, t
    DOUBLE PRECISION :: xlow1, ylow1, zlow1, dx1, dy1, dz1, t1
    DOUBLE PRECISION :: xlow2, ylow2, zlow2, dx2, dy2, dz2, t2

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: q2
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: q1, qc

    CHARACTER(100) :: fname1, fname2, fname3
    CHARACTER(100) :: dir1_fname1, dir2_fname1, dir3_fname1
    CHARACTER(100) :: dir1_fname2, dir2_fname2, dir3_fname2
    character(100) :: dir3_fname3

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: max_grids
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: global_max
    INTEGER, DIMENSION(:), ALLOCATABLE :: grid_max

    !! CHARACTER :: c
    INTEGER :: nstp, ipos, idigit, l1, l2

    INTEGER :: num_args, ix
    CHARACTER(len=12), DIMENSION(:), ALLOCATABLE :: args
    !! Assume file name of the form save_proc<N> (length 10)
    CHARACTER(100) :: dir1, dir2, dir3
    INTEGER :: d1, d2

    logical :: forestclaw1, forestclaw2

    forestclaw1 = .true.
    forestclaw2 = .true.

    num_args = COMMAND_ARGUMENT_COUNT()
    ALLOCATE(args(num_args))  ! I've omitted checking the return status of the allocation

    IF (num_args > 0) THEN
        !! Read from the command line
        DO ix = 1, num_args
            CALL GET_COMMAND_ARGUMENT(ix,args(ix))
        END DO

        READ(args(1),*) d1
        READ(args(2),*) d2
        READ(args(3),*) iframe

        WRITE(dir1,'(A,I0.5)') 'run_',d1
        WRITE(dir2,'(A,I0.5)') 'run_',d2
        WRITE(dir3,'(A,I0.5,A,I0.5)') 'compare_',d1,'_',d2

    ELSE

        !! Size of compare.dat grid
        WRITE(6,'(A)') 'Enter directory 1 (2d) : '
        READ(5,*) dir1

        WRITE(6,'(A)') 'Enter directory 2 (3dx) : '
        READ(5,*) dir2

        WRITE(6,'(A)') 'Enter directory for results (directory must exist) : '
        READ(5,*) dir3

        WRITE(6,'(A)') 'Enter 3d z-level to compare : '
        READ(5,*) zlevel

        WRITE(6,'(A)') 'Enter Frame number to compare : '
        READ(5,*) iframe

    END IF

    !! OPEN(10,file='compare_results.dat')
    !! READ(10,*) iframe
    !! READ(10,*) dir1
    !! READ(10,*) dir2
    !! READ(10,*) dir3  !! results will be stored here in fort.qXXXX and fort.tXXXX
    !! CLOSE(10)

    fname1 = 'fort.qxxxx'
    fname2 = 'fort.txxxx'
    nstp = iframe
    DO ipos = 10, 7, -1
        idigit = MOD(nstp,10)
        fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
        fname2(ipos:ipos) = CHAR(ICHAR('0') + idigit)
        nstp = nstp / 10
    END DO

    fname3 = 'compare_summary.xxxx'
    nstp = iframe
    DO ipos = 20, 17, -1
       idigit = MOD(nstp,10)
       fname3(ipos:ipos) = CHAR(ICHAR('0') + idigit)
       nstp = nstp/10
    END DO

    dir1_fname1 = TRIM(dir1)//'/'//trim(fname1)
    dir2_fname1 = TRIM(dir2)//'/'//trim(fname1)

    dir1_fname2 = TRIM(dir1)//'/'//trim(fname2)
    dir2_fname2 = TRIM(dir2)//'/'//trim(fname2)


    WRITE(6,*) 'Files that will be compared  :'
    WRITE(6,*) dir1_fname1
    WRITE(6,*) dir2_fname1
    WRITE(6,*) ' '
    WRITE(6,*) dir1_fname2
    WRITE(6,*) dir2_fname2
    WRITE(6,*) ' '

    OPEN(10,file=dir1_fname2)
    READ(10,*) t1
    READ(10,*) meqn1
    READ(10,*) ngrids1
    CLOSE(10)

    OPEN(10,file=dir2_fname2)
    READ(10,*) t2
    READ(10,*) meqn2
    READ(10,*) ngrids2
    CLOSE(10)

    IF (meqn1 .NE. meqn2) THEN
        WRITE(6,*) 'meqn1 .ne. meqn2'
        stop
    ENDIF

    IF (t1 .ne. t2) THEN
        WRITE(6,'(A)') 'WARNING : t1 is not exactly equal to t2'
        WRITE(6,'(A10,E30.20)') 't1', t1
        WRITE(6,'(A10,E30.20)') 't2',t2
        WRITE(6,'(A10,E30.16)') 'diff',ABS(t1-t2)
        WRITE(6,*)
        IF (ABS(t1-t2) > 1e-5) THEN
            WRITE(6,*) 't1 and t2 are not close enough to compare files.'
            STOP
        ENDIF
    ENDIF

    IF (ngrids1 .NE. ngrids2) THEN
        WRITE(6,*) 'ngrids1 .ne. ngrids2'
        STOP
    ENDIF
    t = t1
    ngrids = ngrids1
    meqn = meqn1

    WRITE(6,'(A,F24.8)') 'Time of simulation : t = ',t

    dir3_fname3 = TRIM(dir3)//'/'//TRIM(fname3)
    WRITE(6,*) 'A summary of the results will be stored in : '
    WRITE(6,*) dir3_fname3

    ALLOCATE(max_grids(0:ngrids-1,meqn))

    CALL write_tfile(iframe,t,meqn,ngrids,dir3)
    CALL new_qfile(iframe,dir3)

    ALLOCATE(grid_max(meqn),global_max(meqn))
    grid_max = 0
    global_max = 0

    OPEN(10,FILE=dir1_fname1,ERR=100)
    OPEN(20,FILE=dir2_fname1,ERR=200)
    DO n = 0, ngrids-1
        READ(10,*) ngrid1
        READ(10,*) level1
        if (forestclaw1) then
            READ(10,*) blockno1
            READ(10,*) mpirank1
        else
            ngrid1 = ngrid1 - 1
            level1 = level1 - 1
        endif
        !! Read 2d grid
        READ(10,*) mx1
        READ(10,*) my1
        READ(10,*) xlow1
        READ(10,*) ylow1
        READ(10,*) dx1
        READ(10,*) dy1

        !! Read 3dx grid
        READ(20,*) ngrid2
        READ(20,*) level2
        if (forestclaw2) then
            READ(20,*) blockno2
            READ(20,*) mpirank2
        else
            ngrid2 = ngrid2 - 1
            level2 = level2 - 1
        endif
        READ(20,*) mx2
        READ(20,*) my2
        READ(20,*) mz2
        READ(20,*) xlow2
        READ(20,*) ylow2
        READ(20,*) zlow2
        READ(20,*) dx2
        READ(20,*) dy2
        READ(20,*) dz2

        IF (ngrid1 .NE. ngrid2) THEN
            WRITE(6,*) 'ngrid1 .ne. ngrid2; ', ngrid1, ngrid2
            STOP
        ENDIF

        IF (mx1 .NE. mx2) THEN
            WRITE(6,*) 'mx1 .ne. mx2; ',ngrid, mx1, mx2
            STOP
        ENDIF

        IF (xlow1 .NE. xlow2) THEN
            WRITE(6,*) 'xlow1 .ne. xlow2; ', ngrid, xlow1, xlow2
            STOP
        ENDIF

        IF (ylow1 .NE. ylow2) THEN
            WRITE(6,*) 'ylow1 .ne. ylow2; ', ngrid, ylow1, ylow2
            STOP
        ENDIF

        IF (blockno1 .NE. blockno2) THEN
            WRITE(6,*) 'blockno1 .ne. blockno2; ', ngrid, blockno1, blockno2
            STOP
        ENDIF

        IF (zlevel .LT. 0 .or. zlevel .gt. mz2) THEN
            WRITE(6,*) 'zlevel not in [0,mz]; ', zlevel, mz2
            STOP
        ENDIF


        mx = mx1
        my = my1
        mz = mz2
        level = level1
        ngrid = ngrid1
        mpirank = mpirank1
        blockno = blockno1
        !! mpirank = MAX(mpirank1,mpirank2)
        dx = dx1
        dy = dy1
        dz = dz2
        xlow = xlow1
        ylow = ylow1
        zlow = zlow2

        ALLOCATE(q1(mx,my,meqn))
        ALLOCATE(q2(mx,my,mz,meqn))
        ALLOCATE(qc(mx,my,meqn))

        DO m = 1,meqn
            max_grids(ngrid,m) = 0
        END DO
        DO j = 1,my
            DO i = 1,mx
                READ(10,*) (q1(i,j,m),m = 1,meqn)
            end do
        END DO

        DO k = 1,mz
            DO j = 1,my
                DO i = 1,mx
                    READ(20,*) (q2(i,j,k,m),m = 1,meqn)
                END DO
            END DO
        END DO

        do j = 1,my
            do i = 1,mx
                DO m = 1,meqn
                    qc(i,j,m) = q1(i,j,m)-q2(i,j,zlevel,m)
                    IF (ABS(qc(i,j,m)) .GT. global_max(m)) THEN
                        global_max(m) = ABS(qc(i,j,m))  !! Largest difference on all grids
                        grid_max(m) = ngrid        !! max occurs on this grid
                    END IF
                    max_grids(ngrid,m) = MAX(max_grids(ngrid,m),ABS(qc(i,j,m)))
                END DO
            end do 
        end do


        CALL write_qfile2(mx,my,meqn,xlow,ylow,&
                         dx,dy,qc,iframe,ngrid,level,blockno,mpirank,dir3)

        DEALLOCATE(q1, q2, qc)
    END DO
    CLOSE(10)
    CLOSE(20)

    WRITE(6,*) ' '
    DO m = 1,meqn
       WRITE(6,150) global_max(m), m, grid_max(m)
    END DO

    WRITE(6,*) ' '
    OPEN(30,file=dir3_fname3)
    DO m = 1,meqn
        WRITE(30,150) global_max(m), m, grid_max(m)
    END DO
    WRITE(30,*) ' '

    DO n = 0,ngrids-1
        WRITE(30,'(I5,10E16.4)') n, (max_grids(n,m),m=1,meqn)
    END DO
    CLOSE(30)
130 FORMAT(E30.16,'     ',A)
140 FORMAT(I30,   '     ',A)
150 FORMAT(E30.16, '     occurs in field',I5,' on grid ',I5)

    RETURN
100 WRITE(6,*) 'Could not open file ', dir1_fname1
200 WRITE(6,*) 'Could not open file ', dir2_fname1

END PROGRAM compare_files2d_3dx


SUBROUTINE write_tfile(iframe,time,meqn,ngrids,dir3)
    IMPLICIT NONE

    INTEGER :: iframe,meqn,ngrids,maux
    CHARACTER(100) :: dir3, dir3_fname1, dir3_fname2
    CHARACTER(10) :: fname1
    CHARACTER(10) :: fname2
    DOUBLE PRECISION :: time
    INTEGER :: matunit2, nstp,ipos,idigit

    !!  fname1 = 'fort.qxxxx'
    !!  fname2 = 'fort.txxxx'
    !!  matunit2 = 52
    !!  nstp     = iframe
    !!  DO ipos = 10, 7, -1
    !!     idigit = MOD(nstp,10)
    !!     fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
    !!     fname2(ipos:ipos) = CHAR(ICHAR('0') + idigit)
    !!     nstp = nstp / 10
    !!  ENDDO

    WRITE(fname1,'(A,I0.4)') 'fort.q',iframe
    WRITE(fname2,'(A,I0.4)') 'fort.t',iframe

    dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)
    dir3_fname2 = TRIM(dir3)//'/'//trim(fname2)

    maux = 0
    matunit2 = 10
    OPEN(matunit2,file=dir3_fname2)
    WRITE(matunit2,1000) time,meqn,ngrids,maux
    1000 FORMAT(e18.8,'    time', /, &
       i5,'                 meqn'/, &
       i5,'                 ngrids'/, &
       i5,'                 maux'/,/)

    CLOSE(matunit2)

END SUBROUTINE write_tfile

SUBROUTINE new_qfile(iframe,dir3)
    IMPLICIT NONE

    INTEGER :: iframe
    INTEGER :: nstp,ipos,idigit
    CHARACTER(100) :: fname1
    CHARACTER(100) :: dir3, dir3_fname1

    INTEGER fid_com
    COMMON /com_newqfile/ fid_com

    fname1 = 'fort.qxxxx'
    fid_com = 51
    nstp     = iframe
    DO ipos = 10, 7, -1
        idigit = MOD(nstp,10)
        fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
        nstp = nstp / 10
    END DO

    dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)

    OPEN(unit=fid_com,file=dir3_fname1,status='replace')
    CLOSE(unit=fid_com)
END SUBROUTINE new_qfile


SUBROUTINE write_qfile3(mx,my,mz,meqn,xlower,ylower, zlower,&
     dx,dy,dz,q,iframe,patch_num,level,blockno,mpirank,dir3)

    IMPLICIT NONE

    INTEGER :: meqn,mbc,mx,my,mz,mpirank
    INTEGER :: iframe,patch_num, level, blockno
    DOUBLE PRECISION :: xlower, ylower, zlower, dx, dy, dz

    DOUBLE PRECISION :: q(mx,my,mz,meqn)

    CHARACTER(10) :: fname1
    INTEGER :: matunit1
    INTEGER :: nstp,ipos,idigit
    INTEGER i,j,k,mq
    CHARACTER(100) :: dir3_fname1,dir3

    INTEGER fid_com
    COMMON /com_newqfile/ fid_com

    fname1 = 'fort.qxxxx'
    fid_com = 51
    nstp     = iframe
    DO ipos = 10, 7, -1
        idigit = MOD(nstp,10)
        fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
        nstp = nstp / 10
    END DO

    dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)
    OPEN(fid_com,file=dir3_fname1,position='append');

    WRITE(fid_com,1001) patch_num, level, blockno, mpirank, mx, my, mz
1001 FORMAT(i5,'                 grid_number',/, &
       i5,'                 AMR_level',/, &
       i5,'                 block_number',/, &
       i5,'                 mpi_rank',/, &
       i5,'                 mx',/, &
       i5,'                 my',/, &
       i5,'                 mz')


    WRITE(fid_com,1002) xlower,ylower,zlower, dx,dy, dz
1002 FORMAT(e24.16,'    xlow', /, &
       e24.16,'    ylow', /, &
       e24.16,'    zlow', /, &
       e24.16,'    dx', /, &
       e24.16,'    dy', /, &
       e24.16,'    dz',/)

    IF (meqn .GT. 10) THEN
        !!        # Format statement 109 below will not work.
        WRITE(6,'(A,A)') 'Warning (out2.f) : meqn > 10; ', &
                    'change format statement 109.'
        STOP
    ENDIF

    !!      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
    DO k = 1,mz
        DO j = 1,my
            DO i = 1,mx
                DO mq = 1,meqn
                    IF (ABS(q(i,j,k,mq)) .LT. 1d-99) THEN
                        q(i,j,k,mq) = 0.d0
                    ENDIF
                END DO
                WRITE(fid_com,120) (q(i,j,k,mq),mq=1,meqn)
            ENDDO
            WRITE(fid_com,*) ' '
        END DO
        WRITE(fid_com,*) ' '
    END DO
     !! Print the first field assuming it is the MPI rank (i.e. an integer)
120  FORMAT (10E26.16)

    CLOSE(fid_com)

END SUBROUTINE write_qfile3


SUBROUTINE write_qfile2(mx,my,meqn,xlower,ylower, &
     dx,dy,q,iframe,patch_num,level,blockno,mpirank,dir3)

  IMPLICIT NONE

  INTEGER :: meqn,mbc,mx,my,mpirank
  INTEGER :: iframe,patch_num, level, blockno
  DOUBLE PRECISION :: xlower, ylower, dx, dy

  DOUBLE PRECISION q(mx,my,meqn)

  CHARACTER(10) :: fname1
  INTEGER :: matunit1
  INTEGER :: nstp,ipos,idigit
  INTEGER i,j,mq
  CHARACTER(100) :: dir3_fname1,dir3

  INTEGER fid_com
  COMMON /com_newqfile/ fid_com

  fname1 = 'fort.qxxxx'
  fid_com = 51
  nstp     = iframe
  DO ipos = 10, 7, -1
     idigit = MOD(nstp,10)
     fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     nstp = nstp / 10
  ENDDO

  dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)
  OPEN(fid_com,file=dir3_fname1,position='append');

  WRITE(fid_com,1001) patch_num, level, blockno, mpirank, mx, my
1001 FORMAT(i5,'                 grid_number',/, &
            i5,'                 AMR_level',/, &
            i5,'                 block_number',/, &
            i5,'                 mpi_rank',/, &
            i5,'                 mx',/, &
            i5,'                 my')


  WRITE(fid_com,1002) xlower,ylower,dx,dy
1002 FORMAT(e24.16,'    xlow', /, &
          e24.16,'    ylow', /, &
          e24.16,'    dx', /, &
          e24.16,'    dy',/)

  IF (meqn .GT. 10) THEN
     !!        # Format statement 109 below will not work.
     WRITE(6,'(A,A)') 'Warning (out2.f) : meqn > 10; ', &
          'change format statement 109.'
     STOP
  ENDIF

  !!      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
  DO j = 1,my
     DO i = 1,mx
        DO mq = 1,meqn
           IF (ABS(q(i,j,mq)) .LT. 1d-99) THEN
              q(i,j,mq) = 0.d0
           ENDIF
        END DO
        WRITE(fid_com,120) (q(i,j,mq),mq=1,meqn)
     ENDDO
     WRITE(fid_com,*) ' '
  ENDDO
!! Print the first field assuming it is the MPI rank (i.e. an integer)
120 FORMAT (10E26.16)

  CLOSE(fid_com)

END SUBROUTINE write_qfile2
