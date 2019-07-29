c     # ----------------------------------------------------------
c     # Exchange routines - (i,j,mq) ordering
c     # ----------------------------------------------------------
c     # exchange_face_ghost
c     # exchange_corner_ghost
c     # exchange_phys_corner_ghost
c     # ----------------------------------------------------------


c     # Exchange edge ghost data with neighboring grid at same level.

      subroutine fclaw2d_clawpatch46_fort_copy_face(mx,my,mbc,
     &      meqn,qthis,
     &      qneighbor,
     &      iface,
     &      transform_ptr)

      implicit none

      integer mx,my,mbc,meqn,iface, ftransform(9)
      integer*8 transform_ptr
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq, idir
      integer i1,j1, i2, j2

      idir = iface/2

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      do mq = 1,meqn
         if (idir .eq. 0) then
            do j = 1,my
               do ibc = 1,mbc
c                 # Exchange at low side of 'this' grid in
c                 # x-direction (idir == 0)
                  if (iface .eq. 0) then
                     i1 = 1-ibc
                     j1 = j
                  elseif (iface .eq. 1) then
                     i1 = mx+ibc
                     j1 = j
                  endif
                  call fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,
     &                  transform_ptr)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         else
            do jbc = 1,mbc
               do i = 1,mx
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  if (iface .eq. 2) then
                     i1 = i
                     j1 = 1-jbc
                  elseif (iface .eq. 3) then
                     i1 = i
                     j1 = my+jbc
                  endif
                  call fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,
     &                  transform_ptr)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         endif
      enddo

      end

      subroutine fclaw2d_clawpatch46_fort_copy_corner(mx,my,mbc,meqn,
     &      qthis, qneighbor, this_icorner,transform_ptr)
      implicit none

      integer mx, my, mbc, meqn, this_icorner
      integer*8 transform_ptr
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc
      integer i1, j1, i2, j2

c     # Do exchanges for all corners
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (this_icorner .eq. 0) then
                  i1 = 1-ibc
                  j1 = 1-jbc
               elseif (this_icorner .eq. 1) then
                  i1 = mx+ibc
                  j1 = 1-jbc
               elseif (this_icorner .eq. 2) then
                  i1 = 1 -ibc
                  j1 = my+jbc
               else
                  i1 = mx+ibc
                  j1 = my+jbc
               endif

c              # this routine is not yet complete, but the complete one
c              # can now be dropped in.
               call fclaw2d_clawpatch_transform_corner(i1,j1,i2,j2,
     &               transform_ptr)
               qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
            enddo
         enddo
      enddo
      end
