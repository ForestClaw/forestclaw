      subroutine pwconst_write_tfile(iframe,time,mfields,ngrids,maux)
      implicit none

      integer iframe,mfields,ngrids,maux

      character*10 matname1
      character*10 matname2
      double precision time
      integer matunit1, matunit2, nstp,ipos,idigit

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
      write(matunit2,1000) time,mfields,ngrids,maux
 1000 format(e18.8,'    time', /,
     &      i5,'                 mfields'/,
     &      i5,'                 ngrids'/,
     &      i5,'                 maux'/,/)

      close(matunit2)

      open(unit=matunit1,file=matname1,status='replace')
      close(matunit1)

      end

      subroutine pwconst_write_qfile(meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,iframe,patch_num,level,blockno,
     &      mpirank)

      implicit none

      integer meqn,mbc,mx,my, mpirank
      integer iframe,patch_num, level, blockno
      double precision xlower, ylower,dx,dy

      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      character*10 matname1
      integer matunit1
      integer nstp,ipos,idigit
      integer i,j,mq

      matname1 = 'fort.qxxxx'
      matunit1 = 10
      nstp     = iframe
      do ipos = 10, 7, -1
         idigit = mod(nstp,10)
         matname1(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
      enddo

      open(matunit1,file=matname1,access='append');

      write(matunit1,1001) patch_num,level,blockno,mpirank,mx, my
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

      if (meqn .gt. 5) then
c        # Format statement 109 below will not work.
         write(6,'(A,A)') 'Warning (out2.f) : meqn > 5; ',
     &         'change format statement 109.'
         stop
      endif

c      write(6,*) 'WARNING : (claw_out2.f ) Setting q to 0'
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
  120 format (50E26.16)

      close(matunit1)

      end
