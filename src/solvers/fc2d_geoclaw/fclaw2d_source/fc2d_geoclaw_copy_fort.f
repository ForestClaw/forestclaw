c     # ----------------------------------------------------------
c     # Exchange routines - (mq,i,j) ordering
c     # ----------------------------------------------------------
c     # exchange_face_ghost
c     # exchange_corner_ghost
c     # exchange_phys_corner_ghost
c     # ----------------------------------------------------------


c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine fc2d_geoclaw_fort_copy_face(mx,my,mbc,meqn,qthis,
     &      qneighbor,iface,transform_ptr)
      implicit none

      integer mx,my,mbc,meqn,iface, ftransform(9)
      integer*8 transform_ptr
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)

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
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)

               enddo
            enddo
         endif
      enddo
      end


      subroutine fc2d_geoclaw_fort_copy_corner(mx,my,mbc,meqn,
     &      qthis, qneighbor, this_icorner,transform_ptr)
      implicit none

      integer mx, my, mbc, meqn, this_icorner
      integer*8 transform_ptr
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
               qthis(mq,i1,j1) = qneighbor(mq,i2,j2)
            enddo
         enddo
      enddo
      end


      subroutine fc2d_geoclaw_fort_copy_phys_corner(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iface)
      implicit none

      integer mx, my, mbc, meqn, iface, icorner
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ibc, jbc, mq

c     # Fill in corner ghost cells that overlap the physical boundary. In this
c     case, the corner ghost cells are copied from a face neighbor.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (iface .eq. 1) then
                  if (icorner .eq. 1) then
                     qthis(mq,mx+ibc,jbc-mbc) =
     &                     qneighbor(mq,ibc,jbc-mbc)
                     qneighbor(mq,ibc-mbc,jbc-mbc) =
     &                     qthis(mq,mx+ibc-mbc,jbc-mbc)
                  elseif (icorner .eq. 3) then
                     qthis(mq,mx+ibc,my+jbc) =
     &                     qneighbor(mq,ibc,my+jbc)
                     qneighbor(mq,ibc-mbc,my+jbc) =
     &                     qthis(mq,mx+ibc-mbc,my+jbc)
                  endif
               elseif (iface .eq. 3) then
                  if (icorner .eq. 2) then
                     qthis(mq,ibc-mbc,my+jbc) =
     &                     qneighbor(mq,ibc-mbc,jbc)
                     qneighbor(mq,ibc-mbc,jbc-mbc) =
     &                     qthis(mq,ibc-mbc,my+jbc-mbc)
                  elseif(icorner .eq. 3) then
                     qthis(mq,mx+ibc,my+jbc) =
     &                     qneighbor(mq,mx+ibc,jbc)
                     qneighbor(mq,mx+ibc,jbc-mbc) =
     &                     qthis(mq,mx+ibc,my+jbc-mbc)
                  endif
               endif
            enddo
         enddo
      enddo
      end
