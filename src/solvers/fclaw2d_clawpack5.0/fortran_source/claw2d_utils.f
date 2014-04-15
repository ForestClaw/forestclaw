c    # ------------------------------------------------------------------
c    # Internal boundary conditions
c    # ------------------------------------------------------------------

c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine exchange_face_ghost(mx,my,mbc,meqn,qthis,qneighbor,
     &      iface)
      implicit none

      integer mx,my,mbc,meqn,iface
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn, 1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j,ibc,jbc,mq, idir

      idir = iface/2

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      if (idir .eq. 0) then
         do j = 1,my
            do ibc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # x-direction (idir == 0)
                  if (iface .eq. 0) then
                     qthis(mq,1-ibc,j) = qneighbor(mq, mx+1-ibc,j)
                  elseif (iface .eq. 1) then
                     qthis(mq,mx+ibc,j) = qneighbor(mq,ibc,j)
                  endif
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
                     qthis(mq,i,1-jbc) = qneighbor(mq,i,my+1-jbc)
                  elseif (iface .eq. 3) then
                     qthis(mq,i,my+jbc) = qneighbor(mq,i,jbc)
                  endif
               enddo
            enddo
         enddo
      endif
      end


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iside'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost(mx,my,mbc,meqn,qcoarse,
     &      qfine,idir,iface_coarse,num_neighbors,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision sum

      integer mq,r2
      integer i, ic_add, ibc, ii, ifine
      integer j, jc_add, jbc, jj, jfine

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
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (j-1)*refratio + jj
                        if (iface_coarse .eq. 0) then
                           sum = sum + qfine(mq,mx-ifine+1,jfine)
                        elseif (iface_coarse .eq. 1) then
                           sum = sum + qfine(mq,ifine,jfine)
                        endif
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     qcoarse(mq,1-ibc,j+jc_add) = sum/r2
                  else
                     qcoarse(mq,mx+ibc,j+jc_add) = sum/r2
                  endif
               enddo
            enddo
         else
            ic_add = igrid*mx/num_neighbors
            do i = 1,mx/num_neighbors
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (i-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (iface_coarse .eq. 2) then
                           sum = sum + qfine(mq,ifine,my-jfine+1)
                        else
                           sum = sum + qfine(mq,ifine,jfine)
                        endif
                     enddo
                  enddo
                  if (iface_coarse .eq. 2) then
                     qcoarse(mq,i+ic_add,1-jbc) = sum/r2
                  else
                     qcoarse(mq,i+ic_add,my+jbc) = sum/r2
                  endif
               enddo
            enddo
         endif
      enddo

      end

      subroutine interpolate_face_ghost(mx,my,mbc,meqn,qcoarse,qfine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
            j2 = my/num_neighbors + (1-igrid)

            jc_add = igrid*my/num_neighbors
            if (iface_coarse .eq. 0) then
               ic = 1
            elseif (iface_coarse .eq. 1) then
               ic = mx
            endif
            do j = j1, j2
               jc = j + jc_add
               qc = qcoarse(mq,ic,jc)
c              # Compute limited slopes in both x and y. Note we are not
c              # really computing slopes, but rather just differences.
c              # Scaling is accounted for in 'shiftx' and 'shifty', below.
               sl = (qc - qcoarse(mq,ic-1,jc))
               sr = (qcoarse(mq,ic+1,jc) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(mq,ic,jc-1))
               sr = (qcoarse(mq,ic,jc+1) - qc)
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
                        qfine(mq,mx+ifine,jfine) = value
                     elseif (iface_coarse .eq. 1) then
c                       # qfine is at right edge of coarse grid
                        qfine(mq,1-ifine,jfine) = value
                     endif
                  enddo
               enddo
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
               qc = qcoarse(mq,ic,jc)

               sl = (qc - qcoarse(mq,ic-1,jc))
               sr = (qcoarse(mq,ic+1,jc) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(mq,ic,jc-1))
               sr = (qcoarse(mq,ic,jc+1) - qc)
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
                        qfine(mq,ifine,my+jfine) = value
                     else if (iface_coarse .eq. 3) then
c                       # qfine is at top edge of coarse grid
                        qfine(mq,ifine,1-jfine) = value
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
      end

      subroutine exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner_this)
      implicit none

      integer mx, my, mbc, meqn, icorner_this
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq, ibc, jbc

c     # Do exchanges for all corners
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Exchange initiated only at high side (1,3) corners
               if (icorner_this .eq. 0) then
                  qthis(mq,1-ibc,1-jbc) =
     &                  qneighbor(mq,mx+1-ibc,my+1-jbc)
               else if (icorner_this .eq. 1) then
                  qthis(mq,mx+ibc,1-jbc) =
     &                  qneighbor(mq,ibc,my+1-jbc)
               elseif (icorner_this .eq. 2) then
                  qthis(mq,1-ibc,my+jbc) =
     &                  qneighbor(mq,mx+1-ibc,jbc)
               elseif (icorner_this .eq. 3) then
                  qthis(mq,mx+ibc,my+jbc) =
     &                  qneighbor(mq,ibc,jbc)
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
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
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
                        sum = sum + qfine(mq,mx+1-ifine,my+1-jfine)
                     elseif (icorner_coarse .eq. 1) then
                        sum = sum + qfine(mq,ifine,my+1-jfine)
                     elseif (icorner_coarse .eq. 2) then
                        sum = sum + qfine(mq,mx+1-ifine,jfine)
                     elseif (icorner_coarse .eq. 3) then
                        sum = sum + qfine(mq,ifine,jfine)
                     endif
                  enddo
               enddo
               if (icorner_coarse .eq. 0) then
                  qcoarse(mq,1-ibc,1-jbc) = sum/r2
               elseif (icorner_coarse .eq. 1) then
                  qcoarse(mq,mx+ibc,1-jbc) = sum/r2
               elseif (icorner_coarse .eq. 2) then
                  qcoarse(mq,1-ibc,my+jbc) = sum/r2
               elseif (icorner_coarse .eq. 3) then
                  qcoarse(mq,mx+ibc,my+jbc) = sum/r2
               endif
            enddo
         enddo
      enddo


      end

      subroutine interpolate_corner_ghost(mx,my,mbc,meqn,refratio,
     &      qcoarse,qfine,icorner_coarse)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
         qc = qcoarse(mq,ic,jc)

c        # Interpolate coarse grid corners to fine grid corner ghost cells

c        # Compute limited slopes in both x and y. Note we are not
c        # really computing slopes, but rather just differences.
c        # Scaling is accounted for in 'shiftx' and 'shifty', below.
         sl = (qc - qcoarse(mq,ic-1,jc))
         sr = (qcoarse(mq,ic+1,jc) - qc)
         gradx = compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(mq,ic,jc-1))
         sr = (qcoarse(mq,ic,jc+1) - qc)
         grady = compute_slopes(sl,sr,mth)

c        # Loop over fine grid ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = qc + shiftx*gradx + shifty*grady
               if (icorner_coarse .eq. 0) then
                  qfine(mq,mx+ibc,my+jbc) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(mq,1-ibc,my+jbc) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(mq,mx+ibc,1-jbc) = value
               else
                  qfine(mq,1-ibc,1-jbc) = value
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
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer i,j

c     # Do something here....

      end


      subroutine exchange_phys_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iside)
      implicit none

      integer mx, my, mbc, meqn, iside, icorner
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ibc, jbc, mq

c     # Fill in corner ghost cells that overlap the physical boundary. In this
c     case, the corner ghost cells are copied from a face neighbor.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (iside .eq. 1) then
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
               elseif (iside .eq. 3) then
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


c     # Use bc2 to assign boundary conditions to corners with the physical
c     # boundary of another patch.   Use this routine instead of the exchange,
c     above, in the case where q and its neighbor are not at the same level.
c     In this case, we can't just copy, and intepolation/averaging doesn't
c     make sense, since we have a physical boundary condition
      subroutine set_phys_interior_corner_ghost(mx,my,mbc,meqn, q,
     &      icorner, iface, mthbc)
      implicit none

      integer mx, my, mbc, meqn, icorner, iface, mthbc
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
               qc = qcoarse(mq,ic,jc)

c              # Compute limited slopes in both x and y. Note we are not
c              # really computing slopes, but rather just differences.
c              # Scaling is accounted for in 'shiftx' and 'shifty', below.
               sl = (qc - qcoarse(mq,ic-1,jc))
               sr = (qcoarse(mq,ic+1,jc) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(mq,ic,jc-1))
               sr = (qcoarse(mq,ic,jc+1) - qc)
               grady = compute_slopes(sl,sr,mth)

c              # Fill in refined values on coarse grid cell (ic,jc)
               do ii = 1,refratio
                  do jj = 1,refratio
                     shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                     shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                     qfine(mq,(i-1)*refratio + ii,(j-1)*refratio + jj)
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
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
                     sum = sum + qfine(mq,ifine,jfine)
                  enddo
               enddo
               qcoarse(mq,i+ic_add,j + jc_add) = sum/r2
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
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j

      sum = 0
      do i = 1,mx
         do j = 1,my
            sum = sum + q(1,i,j)*dx*dy
         enddo
      enddo
      end
