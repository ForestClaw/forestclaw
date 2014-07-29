c    # ------------------------------------------------------------------
c    # Internal boundary conditions
c    # ------------------------------------------------------------------

c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine exchange_face_ghost(mx,my,mbc,meqn,qthis,
     &      qneighbor,iface,transform_cptr)
      implicit none

      integer mx,my,mbc,meqn,iface, ftransform(9)
      integer*8 transform_cptr
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq, idir
      integer i1,j1, i2(0),j2(0)

      idir = iface/2

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      if (idir .eq. 0) then
         do j = 1,my
            do ibc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at low side of 'this' grid in
c                 # x-direction (idir == 0)
                  if (iface .eq. 0) then
                     i1 = 1-ibc
                     j1 = j
                  elseif (iface .eq. 1) then
                     i1 = mx+mbc
                     j1 = j
                  endif
                  call idxfunc(transform_cptr,i1,j1,i2,j2)
                  qthis(i1,j1,mq) = qneighbor(i2(0),j2(0),mq)

c                  if (iface .eq. 0) then
c                     i2 = ifunc(1-ibc,j,iface,this_block)
c                     j2 = jfunc(1-ibc,j,iface,this_block)
c                     qthis(1-ibc,j,mq) = qneighbor(i2,j2,mq)
c
cc                    # Original code
cc                     qthis(1-ibc,j,mq) = qneighbor(mx+1-ibc,j,mq)
cc                     qneighbor(mx+ibc,j,mq) = qthis(ibc,j,mq)
c                  elseif (iface .eq. 1) then
c    1                i2 = ifunc(mx+ibc,j,iface,this_block)
c    1                j2 = jfunc(mx+ibc,j,iface,this_block)
c                     qthis(mx+ibc,j,mq) = qneighbor(i2,j2,mq)
c
cc                    # Original code
cc                     qthis(mx+ibc,j,mq) = qneighbor(ibc,j,mq)
cc                     qneighbor(1-ibc,j,mq) = qthis(mx+1-ibc,j,mq)
c                  endif
               enddo
            enddo
         enddo
      else
         do i = 1,mx
            do jbc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  if (iface .eq. 2) then
                     i1 = i
                     j1 = 1-jbc
                  elseif (iface .eq. 3) then
                     i1 = i
                     j1 = my+jbc
                  endif
                  call idxfunc(transform_cptr,i1,j1,i2,j2)
                  qthis(i1,j1,mq) = qneighbor(i2(0),j2(0),mq)

c                  if (iface .eq. 2) then
c                     i2 = ifunc(i,1-jbc,iface,this_block)
c                     j2 = jfunc(i,1-jbc,iface,this_block)
c                     qthis(i,1-jbc,mq) = qneighbor(i2,j2,mq)
c
cc                    # Original code
cc                     qthis(i,1-jbc,mq) = qneighbor(i,my+1-jbc,mq)
cc                     qneighbor(i,my+jbc,mq) = qthis(i,jbc,mq)
c                  elseif (iface .eq. 3) then
c                     i2 = ifunc(i,my+jbc,iface,this_block)
c                     j2 = jfunc(i,my+jbc,iface,this_block)
c                     qthis(i,1-jbc,mq) = qneighbor(i2,j2,mq)
c
cc                    # Original code
cc                     qthis(i,my+jbc,mq) = qneighbor(i,jbc,mq)
cc                     qneighbor(i,1-jbc,mq) = qthis(i,my+1-jbc,mq)
c                  endif
               enddo
            enddo
         enddo
      endif
      end


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iface_coarse'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost(mx,my,mbc,meqn,qcoarse,
     &      qfine,idir,iface_coarse,num_neighbors,refratio,igrid,
     &      transform_cptr)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer*8 transform_cptr
      integer num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer mq,r2, m
      integer i, ic_add, ibc, ii, ifine
      integer j, jc_add, jbc, jj, jfine
      integer i1,j1,i2(0:refratio*refratio),j2(0:refratio*refratio)

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
      do mq = 1,meqn
         if (idir .eq. 0) then
            jc_add = igrid*my/num_neighbors
            do j = 1,my/num_neighbors
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer

c                 # Original code
c                 sum = 0
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        ifine = (ibc-1)*refratio + ii
c                        jfine = (j-1)*refratio + jj
c                        if (iface_coarse .eq. 0) then
c                           sum = sum + qfine(mx-ifine+1,jfine,mq)
c                        elseif (iface_coarse .eq. 1) then
c                           sum = sum + qfine(ifine,jfine,mq)
c                        endif
c                     enddo
c                  enddo

                  if (iface_coarse .eq. 0) then
                     i1 = 1-ibc
                     i2 = j+jc_add
                  elseif (iface_coarse .eq. 1) then
                     i1 = mx+ibc
                     i2 = j+jc_add
                  endif
                  call idxfunc(transform_cptr,i1,i2,i2,j2)
                  sum = 0
                  do m = 0,r2-1
                     sum = sum + qfine(i2(m),j2(m),mq)
                  enddo
                  qcoarse(i1,j1,mq) = sum/r2
               enddo
            enddo
         else
            ic_add = igrid*mx/num_neighbors
            do i = 1,mx/num_neighbors
               do jbc = 1,mbc

                  if (iface_coarse .eq. 2) then
                     i1 = i+ic_add
                     i2 = 1-jbc
                  elseif (iface_coarse .eq. 3) then
                     i1 = i+ic_add
                     i2 = mx+jbc
                  endif
                  call idxfunc(transform_cptr,i1,j1,i2,j2)
                  sum = 0
                  do m = 0,3
                     sum = sum + qfine(i2(m),j2(m),mq)
                  enddo
                  qcoarse(i1,j1,mq) = sum/r2

c                  sum = 0
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        ifine = (i-1)*refratio + ii
c                        jfine = (jbc-1)*refratio + jj
c                        if (iface_coarse .eq. 2) then
c                           sum = sum + qfine(ifine,my-jfine+1,mq)
c                        else
c                           sum = sum + qfine(ifine,jfine,mq)
c                        endif
c                     enddo
c                  enddo
c                  if (iface_coarse .eq. 2) then
c                     qcoarse(i+ic_add,1-jbc,mq) = sum/r2
c                  else
c                     qcoarse(i+ic_add,my+jbc,mq) = sum/r2
c                  endif
               enddo
            enddo
         endif
      enddo

      end

      subroutine interpolate_face_ghost(mx,my,mbc,meqn,qcoarse,qfine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid,
     &      transform_cptr)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      integer*8 transform_cptr
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2, m
      integer i, ic1, ic2, ibc, ii, ifine, i1, i2(0:refratio*refratio-1)
      integer j, jc1, jc2, jbc, jj, jfine, j1, j2(0:refratio*refratio-1)
      integer ic_add, jc_add, ic, jc, mth
      double precision shiftx, shifty, gradx, grady, qc, sl, sr, value
      double precision compute_slopes

      mth = 5
      r2 = refratio*refratio

      do mq = 1,meqn
         if (idir .eq. 0) then
c           # this ensures that we get 'hanging' corners

            jc1 = 1-igrid
            jc2 = my/num_neighbors + (1-igrid)

            jc_add = igrid*my/num_neighbors
            if (iface_coarse .eq. 0) then
               ic = 1
            elseif (iface_coarse .eq. 1) then
               ic = mx
            endif
            do j = jc1, jc2
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

               call idxfunc(transform_cptr,ic,jc,i2,j2)
               do m = 0,r2-1
                  shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                  shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                  value = qc + shiftx*gradx + shifty*grady
                  qfine(i2(m),j2(m),mq) = value
               enddo


c               do ibc = 1,mbc
c                  do jj = 1,refratio
cc                    # Fill in interpolated values on fine grid cell
c                     shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
c                     shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
c
c                     value = qc + shiftx*gradx + shifty*grady
c                     i1 = ic
c                     j1 = jc
c                     qfine(i2,j2,mq) = value
c
cc                     ifine = ibc
cc                     jfine = (j-1)*refratio + jj
cc                     if (iface_coarse .eq. 0) then
ccc                       # qfine is at left edge of coarse grid
cc                        qfine(mx+ifine,jfine,mq) = value
cc                     elseif (iface_coarse .eq. 1) then
ccc                       # qfine is at right edge of coarse grid
cc                        qfine(1-ifine,jfine,mq) = value
cc                     endif
c
c                  enddo
c               enddo
            enddo
         else
            ic_add = igrid*mx/num_neighbors
c           # this ensures that we get 'hanging' corners
            i1 = 1-igrid
            i2 = mx/num_neighbors + (1-igrid)

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

               call idxfunc(transform_cptr,ic,jc,i2,j2)
               do m = 0,r2-1
                  shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                  shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                  value = qc + shiftx*gradx + shifty*grady
                  qfine(i2(m),j2(m),mq) = value
               enddo


c               do jbc = 1,mbc
c                  do ii = 1,refratio
cc                    # Fill in interpolated values on fine grid cell
c                     shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
c                     shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio
c
c                     value = (qc + shiftx*gradx + shifty*grady)
c
c                     ifine = (i-1)*refratio + ii
c                     jfine = jbc
c                     if (iface_coarse .eq. 2) then
cc                       # qfine is at bottom edge of coarse grid
c                        qfine(ifine,my+jfine,mq) = value
c                     else if (iface_coarse .eq. 3) then
cc                       # qfine is at top edge of coarse grid
c                        qfine(ifine,1-jfine,mq) = value
c                     endif
c                  enddo
c               enddo

            enddo
         endif
      enddo

      end

      subroutine exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner_this)
      implicit none

      integer mx, my, mbc, meqn, icorner_this
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

c     # Do exchanges for all corners
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Exchange initiated only at high side (1,3) corners
               if (icorner_this .eq. 0) then
                  qthis(1-ibc,1-jbc,mq) =
     &                  qneighbor(mx+1-ibc,my+1-jbc,mq)
               else if (icorner_this .eq. 1) then
                  qthis(mx+ibc,1-jbc,mq) =
     &                  qneighbor(ibc,my+1-jbc,mq)
               elseif (icorner_this .eq. 2) then
                  qthis(1-ibc,my+jbc,mq) =
     &                  qneighbor(mx+1-ibc,jbc,mq)
               elseif (icorner_this .eq. 3) then
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(ibc,jbc,mq)
               endif
            enddo
         enddo
      enddo
      end



c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_corner_ghost(mx,my,mbc,meqn,
     &      refratio,qcoarse,qfine,icorner_coarse)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner_coarse
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2
      integer ifine, jfine

      integer get_block_idx, get_patch_idx
      logical debug_is_on


      r2 = refratio*refratio
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Average fine grid corners onto coarse grid ghost corners
               sum = 0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (ibc-1)*refratio + ii
                     jfine = (jbc-1)*refratio + jj
                     if (icorner_coarse .eq. 0) then
                        sum = sum + qfine(mx+1-ifine,my+1-jfine,mq)
                     elseif (icorner_coarse .eq. 1) then
                        sum = sum + qfine(ifine,my+1-jfine,mq)
                     elseif (icorner_coarse .eq. 2) then
                        sum = sum + qfine(mx+1-ifine,jfine,mq)
                     elseif (icorner_coarse .eq. 3) then
                        sum = sum + qfine(ifine,jfine,mq)
                     endif
                  enddo
               enddo
               if (icorner_coarse .eq. 0) then
                  qcoarse(1-ibc,1-jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 1) then
                  qcoarse(mx+ibc,1-jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 2) then
                  qcoarse(1-ibc,my+jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 3) then
                  qcoarse(mx+ibc,my+jbc,mq) = sum/r2
               endif
            enddo
         enddo
      enddo


      end

      subroutine interpolate_corner_ghost(mx,my,mbc,meqn,refratio,
     &      qcoarse,qfine,icorner_coarse)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ic, jc, mq, ibc,jbc, mth,i,j
      double precision qc, sl, sr, gradx, grady, shiftx, shifty
      double precision compute_slopes, value

      integer get_patch_idx, pidx
      logical debug_is_on

      pidx = get_patch_idx()

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
                  qfine(mx+ibc,my+jbc,mq) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(1-ibc,my+jbc,mq) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(mx+ibc,1-jbc,mq) = value
               else
                  qfine(1-ibc,1-jbc,mq) = value
               endif
            enddo
         enddo
      enddo

      end

c     # ----------------------------------------------------------------------
c     # Physical boundary conditions
c     # ----------------------------------------------------------------------
      subroutine set_phys_corner_ghost(mx,my,mbc,meqn,q,icorner,t,dt,
     &      mthbc)
      implicit none

      integer mx,my,mbc,meqn,icorner, mthbc(4)
      double precision t,dt
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

c     # Do something here....

      end


      subroutine exchange_phys_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iside)
      implicit none

      integer mx, my, mbc, meqn, iside, icorner
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ibc, jbc, mq

c     # Fill in corner ghost cells that overlap the physical boundary. In this
c     case, the corner ghost cells are copied from a face neighbor.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (iside .eq. 1) then
                  if (icorner .eq. 1) then
                     qthis(mx+ibc,jbc-mbc,mq) =
     &                     qneighbor(ibc,jbc-mbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(mx+ibc-mbc,jbc-mbc,mq)
                  elseif (icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(ibc,my+jbc,mq)
                     qneighbor(ibc-mbc,my+jbc,mq) =
     &                     qthis(mx+ibc-mbc,my+jbc,mq)
                  endif
               elseif (iside .eq. 3) then
                  if (icorner .eq. 2) then
                     qthis(ibc-mbc,my+jbc,mq) =
     &                     qneighbor(ibc-mbc,jbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(ibc-mbc,my+jbc-mbc,mq)
                  elseif(icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(mx+ibc,jbc,mq)
                     qneighbor(mx+ibc,jbc-mbc,mq) =
     &                     qthis(mx+ibc,my+jbc-mbc,mq)
                  endif
               endif
            enddo
         enddo
      enddo
      end


c     # Use bc2 to assign boundary conditions to corners with the physical
c     # boundary of another patch.   Use this routine instead of the exchange,
c     above, in the case where q and its neighbor are not at the same level.
c     In this case, we can't just copy, and intepolation/averaging doesn't
c     make sense, since we have a physical boundary condition
      subroutine set_phys_interior_corner_ghost(mx,my,mbc,meqn, q,
     &      icorner, iface, mthbc)
      implicit none

      integer mx, my, mbc, meqn, icorner, iface, mthbc
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer this_side_bc, mthbc_tmp

c     # Redo this side, knowing that it now has valid interior ghost cell data
c     # from which we can copy.

      end

c    # ---------------------------------------------------------------------
c    # Tagging for refinement/coarsening
c    # ---------------------------------------------------------------------


c     # Conservative intepolation to fine grid patch
      subroutine interpolate_to_fine_patch(mx,my,mbc,meqn, qcoarse,
     &      qfine, p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      integer ic,jc,ic_add, jc_add
      double precision qc, shiftx, shifty, sl, sr, gradx, grady
      double precision compute_slopes


c     # Use limiting done in AMRClaw.
      mth = 5

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

      i1 = 1-ig
      i2 = mx/p4est_refineFactor + (1-ig)
      ic_add = ig*mx/p4est_refineFactor

      j1 = 1-jg
      j2 = my/p4est_refineFactor + (1-jg)
      jc_add = jg*my/p4est_refineFactor

      do mq = 1,meqn
         do i = i1,i2
            do j = j1,j2
               ic = i + ic_add
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

c              # Fill in refined values on coarse grid cell (ic,jc)
               do ii = 1,refratio
                  do jj = 1,refratio
                     shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                     shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                     qfine((i-1)*refratio + ii,(j-1)*refratio + jj,mq)
     &                     = qc + shiftx*gradx + shifty*grady
                  enddo
               enddo
            enddo
         enddo
      enddo

      end

      double precision function compute_slopes(sl,sr,mth)
      implicit none

      double precision sl,sr, s, sc, philim, slim
      integer mth

c     # ------------------------------------------------
c     # Slope limiting done in amrclaw - see filpatch.f
c       dupc = valp10 - valc
c       dumc = valc   - valm10
c       ducc = valp10 - valm10
c       du   = dmin1(dabs(dupc),dabs(dumc))        <--
c       du   = dmin1(2.d0*du,.5d0*dabs(ducc))      <-- Not quite sure I follow
c
c       fu = dmax1(0.d0,dsign(1.d0,dupc*dumc))
c
c       dvpc = valp01 - valc
c       dvmc = valc   - valm01
c       dvcc = valp01 - valm01
c       dv   = dmin1(dabs(dvpc),dabs(dvmc))
c       dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
c       fv = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))
c
c       valint = valc + eta1*du*dsign(1.d0,ducc)*fu
c      .      + eta2*dv*dsign(1.d0,dvcc)*fv
c     # ------------------------------------------------

c     # ------------------------------------------------
c     # To see what Chombo does, look in InterpF.ChF
c     # (in Chombo/lib/src/AMRTools) for routine 'interplimit'
c     # Good luck.
c     # ------------------------------------------------

      if (mth .le. 4) then
c        # Use minmod, superbee, etc.
         slim = philim(sl,sr,mth)
         compute_slopes = slim*sl
      else
c        # Use AMRClaw slopes
         sc = (sl + sr)/2.d0
         compute_slopes = min(abs(sl),abs(sr),abs(sc))*
     &         max(0.d0,sign(1.d0,sl*sr))*sign(1.d0,sc)
      endif


      end

      subroutine average_to_coarse_patch(mx,my,mbc,meqn,qcoarse,qfine,
     &      p4est_refineFactor, refratio, igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      integer mq
      double precision sum, r2

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*mx/p4est_refineFactor

      r2 = refratio*refratio
      do mq = 1,meqn
         do j = 1,my/p4est_refineFactor
            do i = 1,mx/p4est_refineFactor
               sum = 0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     sum = sum + qfine(ifine,jfine,mq)
                  enddo
               enddo
               qcoarse(i+ic_add,j + jc_add,mq) = sum/r2
            enddo
         enddo
      enddo
      end

c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
      subroutine compute_sum(mx,my,mbc,meqn,dx,dy,q,sum)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, sum
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j

      sum = 0
      do i = 1,mx
         do j = 1,my
            sum = sum + q(i,j,1)*dx*dy
         enddo
      enddo
      end
