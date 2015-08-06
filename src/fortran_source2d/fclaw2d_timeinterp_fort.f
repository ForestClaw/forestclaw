      subroutine fclaw2d_timeinterp_fort(mx,my,mbc,meqn,psize,
     &      qcurr,qlast,qinterp,alpha,ierror)
      implicit none

      integer mx,my,mbc,meqn,psize, ierror
      double precision alpha
      double precision   qcurr(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision   qlast(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qinterp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, m,mint,k,kfinal

c     # Number of interior layers to compute.  Since we only
c     # interpolate from time-interpolated levels, we only
c     # need two layers.  If we were averaging, we'd need four.
      mint = 2
      ierror = 0
      k = 1

      do m = 1,meqn

c        # Face 0
         do j = 0,my-mint
            do i = 0,mint
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 2
         do j = 0,mint
            do i = mint+1,mx+1
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my+1
            do i = mx-mint+1,mx+1
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 3
         do j = my-mint+1,my+1
            do i = 0,mx-mint
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

      enddo

      kfinal = k-1;
      if (kfinal .ne. psize) then
         ierror = 2
      endif


      end
