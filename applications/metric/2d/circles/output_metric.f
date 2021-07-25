      subroutine output_metric(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower, dx,dy,q,iframe,patch_num,level,
     &      blockno,xp,yp)

      implicit none

      integer maxmx, maxmy,meqn,mbc,mx,my, maux
      integer iframe,patch_num, level, blockno
      double precision xlower, ylower,dx,dy

      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      double precision qv(2)
      integer m,mfields

      character*10 matname1
      integer matunit1
      integer nstp,ipos,idigit
      integer i,j,mq

c     # xp, yp included here
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

      write(matunit1,1001) patch_num, level, blockno, mx, my
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 block_number',/,
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
      mfields = 2
      do j = 1,my
         do i = 1,mx
            qv(1) = curvature(i,j)
            qv(2) = log10(abs(curvature(i,j)-1))
            do m = 1,mfields
               if (abs(qv(m)) .lt. 1d-99) then
                  qv(m) = 0.d0
               endif
            enddo
            write(matunit1,120) 1.0, 1.0
         enddo
         write(matunit1,*) ' '
      enddo
  120 format (5E26.16)

      close(matunit1)

      end
