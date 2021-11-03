c> @file
c> Clawpatch 4.6 timeinterp subroutines

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_timeinterp_t
c>
c> Implementation for clawpack 4.6.
c>    
c> @details @copydetails ::clawpatch_fort_timeinterp_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_timeinterp
     &      (mx,my,mbc,meqn,psize,
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
      mint = mbc
      ierror = 0
      k = 1

c     # Time interpolate only to the interior cells.

      do m = 1,meqn
c        # Face 0
         do j = 1,my-mint
            do i = 1,mint
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 2
         do j = 1,mint
            do i = mint+1,mx
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my
            do i = mx-mint+1,mx
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
               k = k + 1
            enddo
         enddo

c        # Face 3
         do j = my-mint+1,my
            do i = 1,mx-mint
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
