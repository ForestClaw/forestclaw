c     # ----------------------------------------------------------
c>      @file
c>      Exchange routines - (i,j,mq) ordering
c     # ----------------------------------------------------------
c     # exchange_face_ghost
c     # exchange_corner_ghost
c     # exchange_phys_corner_ghost
c     # ----------------------------------------------------------

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_copy_face_t
c>
c> Implementation for clawpack 4.6
c>
c> @details @copydetails ::clawpatch_fort_copy_face_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_copy_face(mx,my,mbc,
     &      meqn,qthis,
     &      qneighbor,
     &      iface,
     &      transform_ptr)

      implicit none

      integer mx,my,mbc,meqn,iface
      integer*8 transform_ptr
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq, idir
      integer i1,j1, i2, j2

      integer a(2,2), f(2)

      idir = iface/2

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      call fclaw2d_clawpatch_build_transform_same(transform_ptr, 
     &                                             a, f)

      if (iface .eq. 0) then
         do mq = 1,meqn
            do j1 = 1,my
               do i1 = 1-mbc,0
c                 # Lower side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         enddo
      else if (iface .eq. 1) then
         do mq = 1,meqn
            do j1 = 1,my
               do i1 = mx+1,mx+mbc
c                 # Upper side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         enddo
      else if (iface .eq. 2) then
         do mq = 1,meqn
            do j1 = 1-mbc,0
               do i1 = 1,mx
c                 # left side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         enddo
      else if (iface .eq. 3) then
         do mq = 1,meqn
            do j1 = my+1,my+mbc
               do i1 = 1,mx
c                 # right side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)
               enddo
            enddo
         enddo
      endif

      end

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_copy_corner_t
c>
c> Implementation for clawpack 4.6
c>
c> @details @copydetails ::clawpatch_fort_copy_corner_t
c--------------------------------------------------------------------
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
