c     # ----------------------------------------------------------------------
c     # This handles the boundary conditions across blocks
c     # All boundaries in the two-patch case are assumed to
c     # oriented identically.
c     # ----------------------------------------------------------------------


c     # Exchange edge ghost data with neighboring grid at same level.
c     # In the multiblock case, 'this' grid always initiates the exchange
c     # at face 'iface'.
      subroutine mb_exchange_face_ghost(mx,my,mbc,meqn,qthis,
     &      qneighbor, iface,iblock)
      implicit none

      integer mx,my,mbc,meqn,iface, iblock
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq

      if (iblock .eq. 1) then
c        # Only initiate exchange from block 0.  I am not sure if the
c        # saves any work, but seems silly to exchange twice.
         return
      endif

c      write(6,*) 'Warning ------------ Setting all BCs to 0 --------'

c      call set_bcs_to_zero(mx,my,mbc,meqn,qthis, qneighbor)


c      return


      do mq = 1,meqn
         if (iface .eq. 0) then
            do j = 1,my
               do ibc = 1,mbc
                  qthis(1-ibc,j,mq) = qneighbor(ibc,j,mq)
                  qneighbor(1-ibc,j,mq) = qthis(ibc,j,mq)
               enddo
            enddo
         elseif (iface .eq. 1) then
            do j = 1,my
               do ibc = 1,mbc
                  qthis(mx+ibc,j,mq) = qneighbor(mx+1-ibc,j,mq)
                  qneighbor(mx+ibc,j,mq) = qthis(mx+1-ibc,j,mq)
               enddo
            enddo
         elseif (iface .eq. 2) then
            do i = 1,mx
               do jbc = 1,mbc
                  qthis(i,1-jbc,mq) = qneighbor(i,jbc,mq)
                   qneighbor(i,1-jbc,mq) = qthis(i,jbc,mq)
               enddo
            enddo
         elseif (iface .eq. 3) then
            do i = 1,mx
               do jbc = 1,mbc
                  qthis(i,my+jbc,mq) = qneighbor(i,my+1-jbc,mq)
                  qneighbor(i,my+jbc,mq) = qthis(i,my+1-jbc,mq)
               enddo
            enddo
         endif
      enddo
      end

      subroutine set_bcs_to_zero(mx,my,mbc,meqn,
     &      qthis, qneighbor)


      integer mx,my,mbc,meqn
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq


      write(6,*) 'Warning ------------ Setting all BCs to 0 --------'

c      return

      do mq = 1,meqn
         do j = 1,my
            do ibc = 1,mbc
               qthis(1-ibc,j,mq) = 0
               qthis(mx+ibc,j,mq) = 0
               qneighbor(1-ibc,j,mq) = 0
               qneighbor(mx+ibc,j,mq) = 0
            enddo
         enddo

         do i = 1-mbc,mx+mbc
            do jbc = 1,mbc
               qthis(i,1-jbc,mq) = 0
               qthis(i,my+ibc,mq) = 0
               qneighbor(i,1-jbc,mq) = 0
               qneighbor(i,my+jbc,mq) = 0
            enddo
         enddo
      enddo

      return
      end



c     # Exchange ghost cells at interior of a block boundary.
      subroutine mb_exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, is_block_bdry, iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer is_block_bdry(0:3)

      integer mq, ibc, jbc, m

c      call set_bcs_to_zero(mx,my,mbc,meqn,qthis, qneighbor)

c      return

      if (iblock .eq. 1) then
c        # Only initiate exchange from block 0.  I am not sure if the
c        # saves any work, but seems silly to exchange twice.
         return
      endif


      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (is_block_bdry(0) .eq. 1) then
                     qthis(1-ibc,1-jbc,mq) =
     &                     qneighbor(ibc,my+1-jbc,mq)
                     qneighbor(1-ibc,my+jbc,mq) =
     &                     qthis(ibc,jbc,mq)
                  else if (is_block_bdry(2) .eq. 1) then
                     qthis(1-ibc,1-jbc,mq) =
     &                     qneighbor(mx+1-ibc,jbc,mq)
                     qneighbor(mx+ibc,1-jbc,mq) =
     &                     qthis(ibc,jbc,mq)
                  endif
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (is_block_bdry(1) .eq. 1) then
                     qthis(mx+ibc,1-jbc,mq) =
     &                     qneighbor(mx+1-ibc,my+1-jbc,mq)
                     qneighbor(mx+ibc,my+jbc,mq) =
     &                     qthis(mx+1-ibc,jbc,mq)
                  elseif (is_block_bdry(2) .eq. 1) then
                     qthis(mx+ibc,1-jbc,mq) =
     &                     qneighbor(ibc,jbc,mq)
                     qneighbor(1-ibc,1-jbc,mq) =
     &                     qthis(mx+1-ibc,jbc,mq)
                  endif
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (is_block_bdry(0) .eq. 1) then
                     qthis(1-ibc,my+ibc,mq) =
     &                     qneighbor(ibc,jbc,mq)
                     qneighbor(1-ibc,1-jbc,mq) =
     &                     qthis(ibc,my+1-jbc,mq)
                  elseif (is_block_bdry(3) .eq. 1) then
                     qthis(1-ibc,my+jbc,mq) =
     &                     qneighbor(mx+1-ibc,my+1-jbc,mq)
                     qneighbor(mx+ibc,my+jbc,mq) =
     &                     qthis(ibc,my+1-jbc,mq)
                  endif
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (is_block_bdry(1) .eq. 1) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(mx+1-ibc,jbc,mq)
                     qneighbor(mx+ibc,1-jbc,mq) =
     &                     qthis(mx+1-ibc,my+1-jbc,mq)
                  elseif (is_block_bdry(3) .eq. 1) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(ibc,my+1-jbc,mq)
                     qneighbor(1-ibc,my+jbc,mq) =
     &                     qthis(mx+1-ibc,my+1-jbc,mq)
                  endif
               enddo
            enddo
         endif
      enddo
      end


c     # Exchange ghost cells at block corner
      subroutine mb_exchange_block_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(1-ibc,1-jbc,mq) = qneighbor(ibc,jbc,mq)
                  qneighbor(1-ibc,1-jbc,mq) = qthis(ibc,jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mx+ibc,1-jbc,mq) = qneighbor(mx+1-ibc,jbc,mq)
                  qneighbor(mx+ibc,1-jbc,mq) = qthis(mx+1-ibc,jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(1-ibc,my+jbc,mq) = qneighbor(ibc,my+1-jbc,mq)
                  qneighbor(1-ibc,my+jbc,mq) = qthis(ibc,my+1-jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(mx+1-ibc,my+1-jbc,mq)
                  qneighbor(mx+ibc,my+jbc,mq) =
     &                  qthis(mx+1-ibc,my+1-jbc,mq)
               enddo
            enddo
         endif
      enddo
      end


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iside'  in direction 'idir' of 'qcoarse'
c     # Assume this is mapped.
      subroutine mb_average_face_ghost(mx,my,mbc,meqn,qcoarse,
     &      qfine,auxcoarse, auxfine,
     &      idir,iface_coarse,p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer p4est_refineFactor
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,1)

      double precision sum, kf, kc, qf

      integer mq,r2
      integer i, ic_add, ibc, ii, ifine
      integer j, jc_add, jbc, jj, jfine

c     # 'iface' is relative to the coarse grid

c      write(6,*)'Mess with mb_average_face_ghost'
c      call set_bcs_to_zero(mx,my,mbc,meqn,qcoarse, qfine)
c      return



      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
      do mq = 1,meqn
         if (idir .eq. 0) then
            jc_add = igrid*my/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (j-1)*refratio + jj
                        if (iface_coarse .eq. 0) then
                           qf = qfine(ifine,jfine,mq)
                           kf = auxfine(ifine,jfine,1)
                        else
                           qf = qfine(mx+1-ifine,jfine,mq)
                           kf = auxfine(mx+1-ifine,jfine,1)
                        endif
                        sum = sum + qf*kf
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     kc = r2*auxcoarse(1-ibc,j+jc_add,1)
                     qcoarse(1-ibc,j+jc_add,mq) = sum/kc
                  else
                     kc = r2*auxcoarse(mx+ibc,j+jc_add,1)
                     qcoarse(mx+ibc,j+jc_add,mq) = sum/kc
                  endif
               enddo
            enddo
         else
            ic_add = igrid*mx/p4est_refineFactor
            do i = 1,mx/p4est_refineFactor
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (i-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (iface_coarse .eq. 2) then
                           qf = qfine(ifine,jfine,mq)
                           kf = auxfine(ifine,jfine,1)
                        else
                           qf = qfine(ifine,my+1-jfine,mq)
                           kf = auxfine(ifine,my+1-jfine,1)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  if (iface_coarse .eq. 2) then
                     kc = r2*auxcoarse(i+ic_add,1-jbc,1)
                     qcoarse(i+ic_add,1-jbc,mq) = sum/kc
                  else
                     kc = r2*auxcoarse(i+ic_add,my+jbc,1)
                     qcoarse(i+ic_add,my+jbc,mq) = sum/kc
                  endif
               enddo
            enddo
         endif
      enddo

      end

      subroutine mb_interpolate_face_ghost(mx,my,mbc,meqn,qcoarse,qfine,
     &      idir,iface_coarse,p4est_refineFactor,refratio,igrid)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer p4est_refineFactor
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2
      integer i, i1, i2, ibc, ii, ifine
      integer j, j1, j2, jbc, jj, jfine
      integer ic_add, jc_add, ic, jc, mth
      double precision shiftx, shifty, gradx, grady, qc, sl, sr, value
      double precision compute_slopes

c     # To be figured out later
      mth = 5

c     # 'iface_coarse is relative to the coarse grid

      do mq = 1,meqn
         if (idir .eq. 0) then
c           # this ensures that we get 'hanging' corners
            j1 = 1-igrid
            j2 = my/p4est_refineFactor + (1-igrid)

            jc_add = igrid*my/p4est_refineFactor
            if (iface_coarse .eq. 0) then
               ic = 1
            elseif (iface_coarse .eq. 1) then
               ic = mx
            endif
            do j = j1, j2
               jc = j + jc_add
               qc = qcoarse(ic,jc,mq)
c              # Compute limited slopes in both x and y. Note we are not
c              # really computing slopes, but rather just differences.
c              # Scaling is accounted for in 'shiftx' and 'shifty', below.
               sl = (qc - qcoarse(ic-1,jc,mq))
               sr = (qcoarse(ic+1,jc,mq) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(ic,jc-1,mq))
               sr = (qcoarse(ic,jc+1,mq) - qc)
               grady = compute_slopes(sl,sr,mth)

               do ibc = 1,mbc
                  do jj = 1,refratio
c                    # Fill in interpolated values on fine grid cell
                     shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
                     shifty = (jj - refratio/2.d0 - 0.5d0)/refratio

                     value = qc + shiftx*gradx + shifty*grady

                     ifine = ibc
                     jfine = (j-1)*refratio + jj
                     if (iface_coarse .eq. 0) then
c                       # qfine is at left edge of coarse grid
                        qfine(1-ifine,jfine,mq) = value
                     elseif (iface_coarse .eq. 1) then
c                       # qfine is at right edge of coarse grid
                        qfine(mx+ifine,jfine,mq) = value
                     endif
                  enddo
               enddo
            enddo
         else
            ic_add = igrid*mx/p4est_refineFactor
c           # this ensures that we get 'hanging' corners
            i1 = 1-igrid
            i2 = mx/p4est_refineFactor + (1-igrid)

            if (iface_coarse .eq. 2) then
               jc = 1
            elseif (iface_coarse .eq. 3) then
               jc = my
            endif
            do i = i1, i2
               ic = i + ic_add
               qc = qcoarse(ic,jc,mq)

               sl = (qc - qcoarse(ic-1,jc,mq))
               sr = (qcoarse(ic+1,jc,mq) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(ic,jc-1,mq))
               sr = (qcoarse(ic,jc+1,mq) - qc)
               grady = compute_slopes(sl,sr,mth)


               do jbc = 1,mbc
                  do ii = 1,refratio
c                    # Fill in interpolated values on fine grid cell
                     shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                     shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

                     value = (qc + shiftx*gradx + shifty*grady)

                     ifine = (i-1)*refratio + ii
                     jfine = jbc
                     if (iface_coarse .eq. 2) then
c                       # qfine is at bottom edge of coarse grid
                        qfine(ifine,1-jfine,mq) = value
                     else if (iface_coarse .eq. 3) then
c                       # qfine is at top edge of coarse grid
                        qfine(ifine,my+jfine,mq) = value
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
      end

c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine mb_average_corner_ghost(mx,my,mbc,meqn,
     &      refratio,qcoarse,qfine,auxcoarse, auxfine,
     &      icorner, is_block_bdry)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner
      integer is_block_bdry(0:3)
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2
      integer ifine, jfine
      double precision kf, qf, kc

      integer get_block_idx, get_patch_idx
      logical debug_is_on

      include "debug.i"

      block_idx = get_block_idx()
      patch_idx = get_patch_idx()
      if (block_idx .eq. 0 .and. patch_idx .eq. 11) then
         call set_debug_on()
      endif

      if (debug_is_on()) then
         write(6,'(A,I2)') 'Before corner averaging (icorner = ',
     &         icorner,')'
         do i = 1-mbc,2
            do j = 1-mbc,2
               write(6,100) i,j,qcoarse(i,j,1)
            enddo
            write(6,*) ' '
         enddo
         write(6,*) ' '
      endif
  100 format(2I5,F24.16)

      r2 = refratio*refratio
      do mq = 1,meqn
         if (icorner .eq. 0) then
            if (debug_is_on()) then
               write(6,*) 'branch : icorner == 0'
            endif
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Average fine grid corners onto coarse grid ghost corners
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (is_block_bdry(0) .eq. 1) then
                           kf = auxfine(ifine,my+1-jfine,1)
                           qf =   qfine(ifine,my+1-jfine,mq)
                        elseif (is_block_bdry(2) .eq. 1) then
                           kf = auxfine(mx+1-ifine,jfine,1)
                           qf =   qfine(mx+1-ifine,jfine,mq)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(1-ibc,1-jbc,1)
                  qcoarse(1-ibc,1-jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Average fine grid corners onto coarse grid ghost corners
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (is_block_bdry(1) .eq. 1) then
                           kf = auxfine(mx+1-ifine,my+1-jfine,1)
                           qf =   qfine(mx+1-ifine,my+1-jfine,mq)
                        elseif (is_block_bdry(2) .eq. 1) then
                           kf = auxfine(ifine,jfine,1)
                           qf =   qfine(ifine,jfine,mq)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(mx+ibc,1-jbc,1)
                  qcoarse(mx+ibc,1-jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Average fine grid corners onto coarse grid ghost corners
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (is_block_bdry(0) .eq. 1) then
                           kf = auxfine(ifine,jfine,1)
                           qf =   qfine(ifine,jfine,mq)
                        elseif (is_block_bdry(3) .eq. 1) then
                           kf = auxfine(mx+1-ifine,my+1-jfine,1)
                           qf =   qfine(mx+1-ifine,my+1-jfine,mq)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(1-ibc,my+jbc,1)
                  qcoarse(1-ibc,my+jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Average fine grid corners onto coarse grid ghost corners
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (is_block_bdry(1) .eq. 1) then
                           kf = auxfine(mx+1-ifine,jfine,1)
                           qf =   qfine(mx+1-ifine,jfine,mq)
                        elseif (is_block_bdry(3) .eq. 1) then
                           kf = auxfine(ifine,my+1-jfine,1)
                           qf =   qfine(ifine,my+1-jfine,mq)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
c                  kc = r2*auxcoarse(mx+ibc,my+jbc,1)
                  kc = 1
                  qcoarse(mx+ibc,my+jbc,mq) = sum/kc
               enddo
            enddo
         endif
      enddo

      if (debug_is_on()) then
         write(6,*) 'After corner averaging'
         do i = 1-mbc,2
            do j = 1-mbc,2
               write(6,100) i,j,qcoarse(i,j,1)
            enddo
            write(6,*) ' '
         enddo
         write(6,*) ' '
      endif

      end

c     # Exchange ghost cells at block corner
      subroutine mb_average_block_corner_ghost(mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, auxcoarse, auxfine,
     &      icorner,iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock, refratio
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,1)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc, ii, jj, r2, ifine, jfine
      double precision sum,qf,kf,kc

c       write(6,*) 'WARNING : (mb_average_block_corner_ghost) : ',
c      &      ' setting to 0'
c
c       do mq = 1,meqn
c          do ibc = 1,mbc
c             do jbc = 1,mbc
c                qcoarse(1-ibc,1-jbc,mq) = 0
c                qcoarse(mx+ibc,my+jbc,mq) = 0
c             enddo
c          enddo
c       enddo
c
c       return


      r2 = refratio**2
      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(ifine,jfine,mq)
                        kf = auxfine(ifine,jfine,1)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(1-ibc,1-ibc,1)
                  qcoarse(1-ibc,1-jbc,mq) = sum/r2
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(mx+1-ifine,jfine,mq)
                        kf = auxfine(mx+1-ifine,jfine,1)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(mx+ibc,1-jbc,mq)
                  qcoarse(mx+ibc,1-jbc,mq) = sum/r2
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(ifine,my+1-jfine,mq)
                        kf = auxfine(ifine,my+1-jfine,mq)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = r2*auxcoarse(1-ibc,my+jbc,1)
                  qcoarse(1-ibc,my+jbc,mq) = sum/r2
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(mx+1-ifine,my+1-jfine,mq)
                        kf = auxfine(mx+1-ifine,my+1-jfine,1)
                        sum = sum + kf*qf
                     enddo
                  enddo
c                  kc = r2*auxcoarse(mx+ibc,my+jbc,1)
                  kc = 1
                  qcoarse(mx+ibc,my+jbc,mq) = sum/r2
               enddo
            enddo
         endif
      enddo
      end


      subroutine mb_interpolate_corner_ghost(mx,my,mbc,meqn,refratio,
     &      qcoarse,qfine,icorner_coarse, is_block_bdry)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer is_block_bdry(0:3)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ic, jc, mq, ibc,jbc, mth,i,j
      double precision qc, sl, sr, gradx, grady, shiftx, shifty
      double precision compute_slopes, value

      mth = 5

      if (icorner_coarse .eq. 0) then
         ic = 1
         jc = 1
      elseif (icorner_coarse .eq. 1) then
         ic = mx
         jc = 1
      elseif (icorner_coarse .eq. 2) then
         ic = 1
         jc = my
      elseif (icorner_coarse .eq. 3) then
         ic = mx
         jc = my
      endif

      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)

c        # Interpolate coarse grid corners to fine grid corner ghost cells

c        # Compute limited slopes in both x and y. Note we are not
c        # really computing slopes, but rather just differences.
c        # Scaling is accounted for in 'shiftx' and 'shifty', below.
         sl = (qc - qcoarse(ic-1,jc,mq))
         sr = (qcoarse(ic+1,jc,mq) - qc)
         gradx = compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(ic,jc-1,mq))
         sr = (qcoarse(ic,jc+1,mq) - qc)
         grady = compute_slopes(sl,sr,mth)

c        # Loop over fine grid ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = qc + shiftx*gradx + shifty*grady
               if (icorner_coarse .eq. 0) then
                  if (is_block_bdry(0) .eq. 1) then
                     qfine(1-ibc,my+jbc,mq) = value
                  elseif (is_block_bdry(2) .eq. 1) then
                     qfine(mx+ibc,1-jbc,mq) = value
                  endif
               elseif (icorner_coarse .eq. 1) then
                  if (is_block_bdry(1) .eq. 1) then
                     qfine(mx+ibc,my+jbc,mq) = value
                  elseif (is_block_bdry(2) .eq. 1) then
                     qfine(1-ibc,1-jbc,mq) = value
                  endif
               elseif (icorner_coarse .eq. 2) then
                  if (is_block_bdry(0) .eq. 1) then
                     qfine(1-ibc,1-jbc,mq) = value
                  elseif (is_block_bdry(3) .eq. 1) then
                     qfine(mx+ibc,my+jbc,mq) = value
                  endif
               elseif (icorner_coarse .eq. 3) then
                  if (is_block_bdry(1) .eq. 1) then
                     qfine(1-ibc,my+jbc,mq) = value
                  elseif (is_block_bdry(3) .eq. 1) then
                     qfine(1-ibc,my+jbc,mq) = value
                  endif
               endif
            enddo
         enddo
      enddo

      end

c     # Exchange ghost cells at block corner
      subroutine mb_interpolate_block_corner_ghost(mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, icorner_coarse, iblock)
      implicit none
      integer mx, my, mbc, meqn, icorner_coarse, iblock, refratio
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc, ii, jj, r2, ifine, jfine
      integer ic, jc, mth
      double precision gradx, grady, shiftx, shifty
      double precision sl, sr, qc, value, compute_slopes

      mth = 5

      if (icorner_coarse .eq. 0) then
         ic = 1
         jc = 1
      elseif (icorner_coarse .eq. 1) then
         ic = mx
         jc = 1
      elseif (icorner_coarse .eq. 2) then
         ic = 1
         jc = my
      elseif (icorner_coarse .eq. 3) then
         ic = mx
         jc = my
      endif

c     # This may not even matter
      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)
         sl = (qc - qcoarse(ic-1,jc,mq))
         sr = (qcoarse(ic+1,jc,mq) - qc)
         gradx = compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(ic,jc-1,mq))
         sr = (qcoarse(ic,jc+1,mq) - qc)
         grady = compute_slopes(sl,sr,mth)

c        # Loop over fine grid ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = qc + shiftx*gradx + shifty*grady
               if (icorner_coarse .eq. 0) then
                  qfine(1-ibc,1-jbc,mq) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(mx+ibc,1-jbc,mq) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(1-ibc,my+jbc,mq) = value
               elseif (icorner_coarse .eq. 3) then
                  qfine(mx+ibc,my+jbc,mq) = value
               endif
            enddo
         enddo
      enddo
      end
