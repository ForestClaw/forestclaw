      subroutine fc3d_clawpack5_fort_timeinterp(mx,my,mz,mbc,meqn,psize,
     &      qcurr,qlast,qinterp,alpha,ierror)
      implicit none

      integer mx,my,mz,mbc,meqn,psize, ierror
      double precision alpha
      double precision   qcurr(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,
     &                         1-mbc:mz+mbc)
      double precision   qlast(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,
     &                         1-mbc:mz+mbc)
      double precision qinterp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,
     &                         1-mbc:mz+mbc)

      integer i,j,k,m,mint,cnt,cntfinal,ng

c     # Number of interior layers to compute.  Since we only
c     # interpolate from time-interpolated levels, we only
c     # need two layers.  If we were averaging, we'd need four.
      mint = 2
      ierror = 0
      cnt = 1

c     # This could also be set to 0 or 1.  But be sure to also set
c     # set ng in fclaw2d_clawpatch_setup_timeinterp.
      ng = 0

      do m = 1,meqn
c        # Face 0
         do k = 1,mz
           do j = 1-ng,my-mint
              do i = 1-ng,mint
                 qinterp(m,i,j,k) = qlast(m,i,j,k) +
     &               alpha*(qcurr(m,i,j,k)-qlast(m,i,j,k))
                 cnt = cnt + 1
              enddo
           enddo

c        # Face 2
           do j = 1-ng,mint
              do i = mint+1,mx+ng
                 qinterp(m,i,j,k) = qlast(m,i,j,k) +
     &               alpha*(qcurr(m,i,j,k)-qlast(m,i,j,k))
                 cnt = cnt + 1
              enddo
           enddo

c        # Face 1
           do j = mint+1,my+ng
              do i = mx-mint+1,mx+ng
                 qinterp(m,i,j,k) = qlast(m,i,j,k) +
     &               alpha*(qcurr(m,i,j,k)-qlast(m,i,j,k))
                 cnt = cnt + 1
              enddo
           enddo

c        # Face 3
           do j = my-mint+1,my+ng
              do i = 1-ng,mx-mint
                 qinterp(m,i,j,k) = qlast(m,i,j,k) +
     &               alpha*(qcurr(m,i,j,k)-qlast(m,i,j,k))
                 cnt = cnt + 1
              enddo
           enddo
        enddo
      enddo

      cntfinal = cnt-1;
      if (cntfinal .ne. psize) then
         ierror = 2
      endif


      end
