c> @file
c> ascii output routines for clawpack 5

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_header_ascii_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::clawpatch_fort_header_ascii_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_header_ascii(matname1,matname2, 
     &     time, meqn,maux, ngrids)
                                                 
      implicit none

      integer meqn,maux,ngrids

      character*11 matname1, matname2
      double precision time
      integer matunit1,matunit2

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
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::clawpatch_fort_output_ascii_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_output_ascii(matname1,
     &      mx,my,meqn,mbc, xlower,ylower, dx,dy,
     &      q,patch_num,level,blockno,mpirank)
      implicit none

      character*11 matname1
      integer meqn,mbc,mx,my
      integer patch_num, level, blockno, mpirank
      double precision xlower, ylower,dx,dy

      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer matunit1
      integer i,j,mq

      matunit1 = 10
      open(matunit1,file=matname1,position='append');

      call fclaw2d_clawpatch5_fort_write_grid_header(matunit1,
     &      mx,my,xlower,ylower, dx,dy,patch_num,level,
     &      blockno,mpirank)


      if (meqn .gt. 20) then
         write(6,'(A,A)') 
     &       'Warning (fclaw2d_clawpatch5_fort_output_ascii) ',
     &         ': meqn > 20;  change format statement 90.'
         stop
      endif

c      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
      do j = 1,my
         do i = 1,mx
            do mq = 1,meqn
               if (abs(q(mq,i,j)) .lt. 1d-99) then
                  q(mq,i,j) = 0.d0
               endif
            enddo
            write(matunit1,120) (q(mq,i,j),mq=1,meqn)
         enddo
         write(matunit1,*) ' '
      enddo
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
      subroutine fclaw2d_clawpatch5_fort_write_grid_header(matunit1,
     &      mx,my,xlower,ylower, dx,dy,patch_num,level,
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
