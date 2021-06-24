      subroutine square5_fort_header_ascii
     &      (matname1,matname2, time,meqn,maux,ngrids)
      implicit none

      integer meqn,ngrids, maux, mfields

      character*11 matname1
      character*11 matname2
      double precision time
      integer matunit1, matunit2

      matunit1 = 10
      matunit2 = 15

      mfields = meqn + 2  !! Include error and exact solution
      open(unit=matunit2,file=matname2)
      write(matunit2,1000) time,mfields,ngrids,maux,2
 1000 format(e30.20,'    time', /,
     &      i5,'                 meqn'/,
     &      i5,'                 ngrids'/,
     &      i5,'                 num_aux'/,
     &      i5,'                 num_dim')

      close(matunit2)

      open(unit=matunit1,file=matname1,status='replace')
      close(matunit1)

      end




      subroutine square5_fort_write_file(matname1,
     &      mx,my,meqn,mbc, xlower,ylower, dx,dy,
     &      q,error,soln, time, patch_num,level,blockno,mpirank)

      implicit none

      character*10 matname1
      integer meqn,mbc,mx,my
      integer patch_num, level, blockno, mpirank
      double precision xlower, ylower,dx,dy,time
      double precision xc,yc,qc

      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision soln(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer matunit1
      integer i,j,mq

c      double precision swirl_divergence, divu

      integer mapping
      common /mapping_comm/ mapping

      integer*8 cont, get_context

      cont = get_context()

      matunit1 = 10
      open(matunit1,file=matname1,position='append');

      call square5_fort_write_grid_header(matunit1,
     &      mx,my,xlower,ylower, dx,dy,patch_num,level,
     &      blockno,mpirank)


      if (meqn .gt. 5) then
         write(6,'(A,A,A)')
     &         'Warning (fclaw2d_fort_write_grid_header.f) ',
     &         ': meqn > 5; change format statement 120.'
         stop
      endif


      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            do mq = 1,meqn
               if (abs(q(mq,i,j)) .lt. 1d-99) then
                  q(mq,i,j) = 0.d0
               endif
            enddo
            qc = soln(1,i,j)
            if (abs(qc) .lt. 1d-99) then
               qc = 0.d0
            endif
            if (abs(error(1,i,j)) .lt. 1d-99) then
               error(1,i,j) = 0.d0
            endif
c            divu = swirl_divergence(xc,yc)
            write(matunit1,120) (q(mq,i,j),mq=1,meqn),qc,
     &            error(1,i,j)
         enddo
         write(matunit1,*) ' '
      enddo
c     # This statement is checked above (meqn <= 5)
  120 format (5E26.16)

      close(matunit1)

      end

      subroutine square5_fort_write_grid_header
     &      (matunit1, mx,my,xlower,ylower, dx,dy,patch_num,level,
     &      blockno,mpirank)

      implicit none

      integer matunit1, mx, my
      integer patch_num, level, blockno, mpirank
      double precision xlower, ylower,dx,dy


      write(matunit1,1001) patch_num, level, blockno, mpirank, mx, my
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 block_number',/,
     &       i5,'                 mpi_rank',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')


      write(matunit1,1002) xlower,ylower,dx,dy
 1002 format(e24.16,'    xlow', /,
     &       e24.16,'    ylow', /,
     &       e24.16,'    dx', /,
     &       e24.16,'    dy',/)


      end

