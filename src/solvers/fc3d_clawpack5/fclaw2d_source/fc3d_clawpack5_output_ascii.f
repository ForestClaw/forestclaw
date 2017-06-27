      subroutine fc3d_clawpack5_fort_write_header(iframe, time, meqn,
     &                                            maux, ngrids)
      implicit none

      integer iframe,meqn,maux,ngrids

      character*11 matname1, matname2
      double precision time
      integer nstp,ipos,idigit,matunit1,matunit2

      matname1 = 'fort.qxxxx'
      matname2 = 'fort.txxxx'
      matunit1 = 10
      matunit2 = 15
      nstp     = iframe
      do ipos = 10, 7, -1
         idigit = mod(nstp,10)
         matname1(ipos:ipos) = char(ichar('0') + idigit)
         matname2(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
      enddo

      open(unit=matunit2,file=matname2)
      write(matunit2,1000) time,meqn,ngrids,maux,2
 1000 format(e30.20,'    time', /,
     &      i5,'                 meqn'/,
     &      i5,'                 ngrids'/,
     &      i5,'                 num_aux'/,
     &      i5,'                 num_dim')

      close(matunit2)

      open(unit=matunit1,file=matname1,status='replace')
      close(matunit1)

      end


      subroutine fc3d_clawpack5_fort_write_file(matname1,
     &      mx,my,mz,meqn,mbc,xlower,ylower,zlower,dx,dy,
     &      dz,q,patch_num,level,blockno,mpirank)
      implicit none

      character*11 matname1
      integer meqn,mbc,mx,my,mz
      integer patch_num, level, blockno, mpirank
      double precision xlower,ylower,zlower,dx,dy,dz

      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

      integer matunit1
      integer i,j,k,mq

      matunit1 = 10
      open(matunit1,file=matname1,position='append');

      call fc3d_clawpack5_fort_write_grid_header(matunit1,
     &      mx,my,mz,xlower,ylower,zlower,dx,dy,dz,patch_num,level,
     &      blockno,mpirank)


      if (meqn .gt. 5) then
         write(6,'(A,A)') 'Warning (fclaw2d_fort_write_file.f) ',
     &         ': meqn > 5;  change format statement 109.'
         stop
      endif

c      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
      do k = 1,mz
        do j = 1,my
           do i = 1,mx
              do mq = 1,meqn
                 if (abs(q(mq,i,j,k)) .lt. 1d-99) then
                    q(mq,i,j,k) = 0.d0
                 endif
              enddo
              write(matunit1,120) (q(mq,i,j,k),mq=1,meqn)
           enddo
           write(matunit1,*) ' '
        enddo
      enddo
  120 format (5E26.16)

      close(matunit1)

      end

      subroutine fc3d_clawpack5_fort_write_grid_header(matunit1,
     &      mx,my,mz,xlower,ylower,zlower,dx,dy,dz,patch_num,level,
     &      blockno,mpirank)

      implicit none

      integer matunit1, mx, my, mz
      integer patch_num, level, blockno, mpirank
      double precision xlower,ylower,zlower,dx,dy,dz


      write(matunit1,1001) patch_num,level,blockno,mpirank,mx,my,mz
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 block_number',/,
     &       i5,'                 mpi_rank',/,
     &       i5,'                 mx',/,
     &       i5,'                 my',/,
     &       i5,'                 mz')


      write(matunit1,1002) xlower,ylower,zlower,dx,dy,dz
 1002 format(e24.16,'    xlow', /,
     &       e24.16,'    ylow', /,
     &       e24.16,'    zlow', /,
     &       e24.16,'    dx', /,
     &       e24.16,'    dy', /,
     &       e24.16,'    dz',/)


      end
