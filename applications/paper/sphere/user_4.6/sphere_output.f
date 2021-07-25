      subroutine sphere_fort_header_ascii
     &      (matname1,matname2, time,meqn,maux,ngrids)
      implicit none

      integer iframe,meqn,ngrids, maux, mfields

      character*11 matname1
      character*11 matname2
      double precision time
      integer matunit1, matunit2

      matunit1 = 10
      matunit2 = 15

      mfields = meqn + 3  !! Include error and exact solution and area
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

      subroutine sphere_fort_write_file(matname1,
     &      mx,my,meqn,mbc, maux, xlower,ylower, dx,dy,
     &      q,error,soln, aux, time, global_patchno,level,
     &      blockno, compute_error, mpirank)
      implicit none

      character*10 matname1
      integer meqn,mbc,mx,my, maux
      integer global_patchno, level, blockno, mpirank
      integer compute_error
      double precision xlower, ylower,dx,dy,time

      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision err, area, qc
      integer matunit1
      integer i,j,mq

      integer mapping
      common /mapping_comm/ mapping

      integer*8 cont, get_context

      cont = get_context()

      matunit1 = 10
      open(matunit1,file=matname1,position='append');

      call swirl46_fort_write_grid_header(matunit1,
     &      mx,my,xlower,ylower, dx,dy,global_patchno,level,
     &      blockno,mpirank)

      if (meqn .gt. 5) then
         write(6,'(A,A,A)')
     &         'Warning (fclaw2d_fort_write_grid_header.f) ',
     &         ': meqn > 5; change format statement 120.'
         stop
      endif

      do j = 1,my
         do i = 1,mx
            do mq = 1,meqn
               if (abs(q(i,j,mq)) .lt. 1d-99) then
                  q(i,j,mq) = 0.d0
               endif
            enddo
            if (compute_error == 1) then 
                qc = soln(i,j,1)
                if (abs(qc) .lt. 1d-99) then
                  qc = 0.d0
                endif
                if (abs(error(i,j,1)) .lt. 1d-99) then
                    error(i,j,1) = 0.d0
                endif
                err = error(i,j,1)
            else
                qc = q(i,j,1)
                err = 0
            end if
            area = aux(i,j,1)*dx*dy
            write(matunit1,120) (q(i,j,mq),mq=1,meqn),qc,
     &            err, area
         enddo
         write(matunit1,*) ' '
      enddo
c     # This statement is checked above (meqn <= 5)
  120 format (5E26.16)

      close(matunit1)

      end

      subroutine swirl46_fort_write_grid_header
     &      (matunit1, mx,my,xlower,ylower, dx,dy,global_patchno,level,
     &      blockno,mpirank)

      implicit none

      integer matunit1, mx, my
      integer global_patchno, level, blockno, mpirank
      double precision xlower, ylower,dx,dy


      write(matunit1,1001) global_patchno, level, blockno, mpirank, 
     &               mx, my
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

