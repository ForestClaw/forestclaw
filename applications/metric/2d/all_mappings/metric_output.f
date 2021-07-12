      subroutine metric_output(meqn,mbc,mx,my,
     &      xlower,ylower, dx,dy,q,iframe,
     &      patch_num,level,blockno,mpirank)
     &      bind(c,name="metric_output")

      implicit none

      integer meqn,mbc,mx,my, maux, mpirank
      integer iframe,patch_num, level, blockno
      double precision xlower, ylower,dx,dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer m,mfields

      character*10 matname1
      integer matunit1
      integer nstp,ipos,idigit
      integer i,j,mq

      include 'fclaw2d_metric_terms.i'

      matname1 = 'fort.qxxxx'
      matunit1 = 10
      nstp     = iframe
      do ipos = 10, 7, -1
         idigit = mod(nstp,10)
         matname1(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
      enddo

      open(matunit1,file=matname1,position='append');

      write(matunit1,1001) patch_num, level, blockno, mpirank,mx, my
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 block_number',/,
     &       i5,'                 mpirank',/,
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
      mfields = 3
      do j = 1,my
         do i = 1,mx
            do m = 1,mfields
               if (abs(q(i,j,m)) .lt. 1d-99) then
                  q(i,j,m) = 0.d0
               endif
            enddo
            write(matunit1,120) (q(i,j,m), m=1,mfields)
         enddo
         write(matunit1,*) ' '
      enddo
  120 format (5E26.16)

      close(matunit1)

      end
