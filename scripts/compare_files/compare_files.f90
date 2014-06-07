!!
!! Compare two fort files in two different directories.
!! To compile :
!!      gfortran -o compare_files compare_files.f90
!!
!! To use : Create a file called

PROGRAM compare_files
  IMPLICIT NONE

  INTEGER :: iframe
  CHARACTER(100) :: dir1, dir2, dir3

  INTEGER :: i,j, n, m

  INTEGER :: ngrids, mfields
  INTEGER :: ngrids1, mfields1
  INTEGER :: ngrids2, mfields2
  INTEGER :: ngrid, level, mx, my, blockno
  INTEGER :: ngrid1, level1, mx1, my1, blockno1
  INTEGER :: ngrid2, level2, mx2, my2, blockno2

  DOUBLE PRECISION :: xlow, ylow, dx, dy, t
  DOUBLE PRECISION :: xlow1, ylow1, dx1, dy1, t1
  DOUBLE PRECISION :: xlow2, ylow2, dx2, dy2, t2

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: q1, q2, qc

  CHARACTER(100) :: fname1, fname2, fname3
  CHARACTER(100) :: dir1_fname1, dir2_fname1, dir3_fname1
  CHARACTER(100) :: dir1_fname2, dir2_fname2, dir3_fname2
  character(100) :: dir3_fname3

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: max_grids
  DOUBLE PRECISION :: global_max
  INTEGER :: grid_max, field_max

  !! CHARACTER :: c
  INTEGER :: nstp, ipos, idigit, l1, l2

!! Size of compare.dat grid
  WRITE(6,'(A)') 'Enter directory 1 : '
  READ(5,*) dir1

  WRITE(6,'(A)') 'Enter directory 2 : '
  READ(5,*) dir2

  WRITE(6,'(A)') 'Enter directory for results (directory must exist) : '
  READ(5,*) dir3

  WRITE(6,'(A)') 'Enter Frame number to compare : '
  READ(5,*) iframe

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
  ENDDO

  fname3 = 'compare_summary.xxxx'
  nstp = iframe
  DO ipos = 20, 17, -1
     idigit = MOD(nstp,10)
     fname3(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     nstp = nstp/10
  ENDDO

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

  dir3_fname3 = TRIM(dir3)//'/'//TRIM(fname3)
  WRITE(6,*) 'A summary of the results will be stored in : '
  WRITE(6,*) dir3_fname3


  OPEN(10,file=dir1_fname2)
  READ(10,*) t1
  READ(10,*) mfields1
  READ(10,*) ngrids1
  CLOSE(10)

  OPEN(10,file=dir2_fname2)
  READ(10,*) t2
  READ(10,*) mfields2
  READ(10,*) ngrids2
  CLOSE(10)

  IF (mfields1 .NE. mfields2) THEN
     WRITE(6,*) 'mfields1 .ne. mfields2'
     stop
  ENDIF

  IF (ABS(t1-t2) > 1e-13) THEN
     WRITE(6,*) 't1 .ne. t2', t1, t2
     stop
  ENDIF

  IF (ngrids1 .NE. ngrids2) THEN
     WRITE(6,*) 'ngrids1 .ne. ngrids2'
     STOP
  ENDIF
  t = t1
  ngrids = ngrids1
  mfields = mfields1

  ALLOCATE(max_grids(0:ngrids-1,mfields))

  CALL write_tfile(iframe,t,mfields,ngrids,dir3)
  CALL new_qfile(iframe,dir3)

  grid_max = 0
  field_max = 0
  global_max = 0

  OPEN(10,file=dir1_fname1)
  OPEN(20,file=dir2_fname1)
  DO n = 0, ngrids-1
     READ(10,*) ngrid1
     READ(10,*) level1
     READ(10,*) blockno1
     READ(10,*) mx1
     READ(10,*) my1
     READ(10,*) xlow1
     READ(10,*) ylow1
     READ(10,*) dx1
     READ(10,*) dy1

     READ(20,*) ngrid2
     READ(20,*) level2
     READ(20,*) blockno2
     READ(20,*) mx2
     READ(20,*) my2
     READ(20,*) xlow2
     READ(20,*) ylow2
     READ(20,*) dx2
     READ(20,*) dy2

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


     mx = mx1
     my = my1
     level = level1
     ngrid = ngrid1
     blockno = blockno1
     dx = dx1
     dy = dy1
     xlow = xlow1
     ylow = ylow1


     ALLOCATE(q1(mx,my,mfields))
     ALLOCATE(q2(mx,my,mfields))
     ALLOCATE(qc(mx,my,mfields))

     DO m = 1,mfields
        max_grids(ngrid,m) = 0
     ENDDO
     DO j = 1,my
        DO i = 1,mx
           READ(10,*) (q1(i,j,m),m = 1,mfields)
           READ(20,*) (q2(i,j,m),m = 1,mfields)
           DO m = 1,mfields
              qc(i,j,m) = ABS(q1(i,j,m)-q2(i,j,m))
              IF (qc(i,j,m) .GT. global_max) THEN
                 global_max = qc(i,j,m)  !! Largest difference on all grids
                 grid_max = ngrid        !! max occurs on this grid
                 field_max = m           !! max occurs in this field
              ENDIF
              max_grids(ngrid,m) = MAX(max_grids(ngrid,m),qc(i,j,m))
           ENDDO
        ENDDO
     ENDDO

     CALL write_qfile(mx,my,mfields,xlow,ylow, &
          dx,dy,qc,iframe,ngrid,level,blockno,dir3)


     DEALLOCATE(q1, q2, qc)
  ENDDO
  CLOSE(10)
  CLOSE(20)

  OPEN(30,file=dir3_fname3)
  WRITE(30,130) global_max, 'Global maximum'
  WRITE(30,130) global_max, 'Global maximum'
  WRITE(30,140) grid_max, 'Global maximum occurs on this grid'
  WRITE(30,140) field_max, 'Global maximum occurs in this field'
  WRITE(30,*) ' '
  DO n = 0,ngrids-1
     WRITE(30,'(I5,10E16.4)') n, (max_grids(ngrid,m),m=1,mfields)
  ENDDO
  CLOSE(30)
130 FORMAT(E30.16,'     ',A)
140 FORMAT(I30,   '     ',A)

END PROGRAM compare_files


SUBROUTINE write_tfile(iframe,time,mfields,ngrids,dir3)
  IMPLICIT NONE

  INTEGER :: iframe,mfields,ngrids,maux
  CHARACTER(100) :: dir3, dir3_fname1, dir3_fname2
  CHARACTER(10) :: fname1
  CHARACTER(10) :: fname2
  DOUBLE PRECISION :: time
  INTEGER :: matunit2, nstp,ipos,idigit

  fname1 = 'fort.qxxxx'
  fname2 = 'fort.txxxx'
  matunit2 = 52
  nstp     = iframe
  DO ipos = 10, 7, -1
     idigit = MOD(nstp,10)
     fname1(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     fname2(ipos:ipos) = CHAR(ICHAR('0') + idigit)
     nstp = nstp / 10
  ENDDO

  dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)
  dir3_fname2 = TRIM(dir3)//'/'//trim(fname2)

  maux = 0
  OPEN(matunit2,file=dir3_fname2)
  WRITE(matunit2,1000) time,mfields,ngrids,maux
1000 FORMAT(e18.8,'    time', /, &
          i5,'                 mfields'/, &
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
  ENDDO

  dir3_fname1 = TRIM(dir3)//'/'//trim(fname1)

  OPEN(unit=fid_com,file=dir3_fname1,status='replace')
  CLOSE(unit=fid_com)
END SUBROUTINE new_qfile


SUBROUTINE write_qfile(mx,my,mfields,xlower,ylower, &
     dx,dy,q,iframe,patch_num,level,blockno,dir3)

  IMPLICIT NONE

  INTEGER :: mfields,mbc,mx,my
  INTEGER :: iframe,patch_num, level, blockno
  DOUBLE PRECISION :: xlower, ylower,dx,dy

  DOUBLE PRECISION q(mx,my,mfields)

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

  WRITE(fid_com,1001) patch_num, level, blockno, mx, my
1001 FORMAT(i5,'                 grid_number',/, &
          i5,'                 AMR_level',/, &
          i5,'                 block_number',/, &
          i5,'                 mx',/, &
          i5,'                 my')


  WRITE(fid_com,1002) xlower,ylower,dx,dy
1002 FORMAT(e24.16,'    xlow', /, &
          e24.16,'    ylow', /, &
          e24.16,'    dx', /, &
          e24.16,'    dy',/)

  IF (mfields .GT. 10) THEN
     !!        # Format statement 109 below will not work.
     WRITE(6,'(A,A)') 'Warning (out2.f) : mfields > 5; ', &
          'change format statement 109.'
     STOP
  ENDIF

  !!      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
  DO j = 1,my
     DO i = 1,mx
        DO mq = 1,mfields
           IF (ABS(q(i,j,mq)) .LT. 1d-99) THEN
              q(i,j,mq) = 0.d0
           ENDIF
        END DO
        WRITE(fid_com,120) INT(q(i,j,1)), (q(i,j,mq),mq=1,mfields-1)
     ENDDO
     WRITE(fid_com,*) ' '
  ENDDO
!! Print the first field assuming it is the MPI rank (i.e. an integer)
120 FORMAT (I5,10E26.16)

  CLOSE(fid_com)

END SUBROUTINE write_qfile
