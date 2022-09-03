c> @file
c> ascii output routines for clawpack 4.6

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_header_ascii_t
c>
c> Implementation for clawpack 4.6.
c>
c> @details @copydetails ::clawpatch_fort_header_ascii_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_header_ascii
     &      (matname1,matname2, time,meqn,maux,ngrids)
      implicit none

      integer meqn,ngrids, maux

      character*11 matname1
      character*11 matname2
      double precision time
      integer matunit1, matunit2

      matunit1 = 10
      matunit2 = 15

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

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_output_ascii_t
c>
c> Implementation for clawpack 4.6.
c>
c> @details @copydetails ::clawpatch_fort_output_ascii_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_output_ascii(matname1,
     &      mx,my,meqn,mbc, xlower,ylower, dx,dy,
     &      q,patch_num,level,blockno,mpirank)

      implicit none

      character(len=11) matname1
      integer meqn,mbc,mx,my
      integer patch_num
      integer level, blockno, mpirank
      double precision xlower, ylower,dx,dy

      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer matunit1
      integer i,j,mq

      matunit1 = 10
      open(matunit1,file=matname1,position='append');

      call fclaw2d_clawpatch46_fort_write_grid_header(matunit1,
     &      mx,my,xlower,ylower, dx,dy,patch_num,level,
     &      blockno,mpirank)


      if (meqn .gt. 20) then
         write(6,'(A,A,A,I5,A)')     
     &         'Warning (fclaw2d_clawpatch46_output_ascii.f) ',
     &         ': meqn > 20; change format statement 120.', 
     &         '(meqn = ',meqn,')'
         stop
      endif

      do j = 1,my
         do i = 1,mx
            do mq = 1,meqn
               if (abs(q(i,j,mq)) .lt. 1d-99) then
                  q(i,j,mq) = 0.d0
               endif
            enddo
            write(matunit1,120) (q(i,j,mq),mq=1,meqn)
         enddo
         write(matunit1,*) ' '
      enddo
c     # This statement is checked above (meqn <= 20)
  120 format (20E26.16)

      close(matunit1)

      end


c--------------------------------------------------------------------
c> @brief Writes the header for a grid
c>
c> @param [in] matunit1 handle for data file
c> @param [in] mx, my the number of cells in the x and y directions
c> @param [in] xlower, ylower lower left coordinate of patch
c> @param [in] dx, dy the spacings in the x and y directions
c> @param [in] patch_num the patch number
c> @param [in] level the level
c> @param [in] blockno the block number
c> @param [in] mpiran the mpi rank
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_write_grid_header
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
