c     # ----------------------------------------------------------
c>      @file
c> .    Exchange routines - (mq,i,j) ordering
c     # ----------------------------------------------------------
c     # exchange_face_ghost
c     # exchange_corner_ghost
c     # exchange_phys_corner_ghost
c     # ----------------------------------------------------------


c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_copy_face_t
c>
c> Implementation for clawpack 5
c>
c> @details @copydetails ::clawpatch_fort_copy_face_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_copy_face(mx,my,mbc,
     &      meqn,qthis,qneighbor,iface,transform_ptr)
      implicit none

      integer mx,my,mbc,meqn,iface
      integer*8 transform_ptr
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq
      integer i1,j1, i2, j2

      integer a(2,2), f(2)

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      call fclaw2d_clawpatch_build_transform_same(transform_ptr, 
     &                                             a, f)

      if (iface .eq. 0) then
         do j1 = 1,my
            do i1 = 1-mbc,0
               do mq = 1,meqn
c                 # Lower side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)
               enddo
            enddo
         enddo
      else if (iface .eq. 1) then
         do j1 = 1,my
            do i1 = mx+1,mx+mbc
               do mq = 1,meqn
c                 # Upper side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)
               enddo
            enddo
         enddo
      else if (iface .eq. 2) then
         do j1 = 1-mbc,0
            do i1 = 1,mx
               do mq = 1,meqn
c                 # left side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)
               enddo
            enddo
         enddo
      else if (iface .eq. 3) then
         do j1 = my+1,my+mbc
            do i1 = 1,mx
               do mq = 1,meqn
c                 # right side
                  i2 = a(1,1)*i1 + a(1,2)*j1 + f(1)
                  j2 = a(2,1)*i1 + a(2,2)*j1 + f(2)
                  qthis(mq,i1,j1) = qneighbor(mq,i2,j2)
               enddo
            enddo
         enddo
      endif


      end


c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_copy_corner_t
c>
c> Implementation for clawpack 5
c>
c> @details @copydetails ::clawpatch_fort_copy_corner_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_copy_corner(mx,my,mbc,meqn,
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


c--------------------------------------------------------------------
c> @brief Copy the corner ghost cells from a face neighbor
c>
c> @param[in] mx, my the number of cells in the x and y directions
c> @param[in] mbc the number of ghost cells
c> @param[in] meqn the number of equations
c> @param[in,out] q this this solution
c> @param[in,out] qneighbor the neighbor the neighbor solution
c> @param[in] icorner the corner to copy values to
c> @param[in] iface the face that the neighbor is on
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_copy_phys_corner(mx,my,
     &      mbc,meqn, qthis, qneighbor, icorner, iface)
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
