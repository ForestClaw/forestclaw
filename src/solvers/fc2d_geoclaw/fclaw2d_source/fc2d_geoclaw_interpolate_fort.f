c     # ----------------------------------------------------------
c     # Interpolation routines - (i,j,mq) ordering
c     # ----------------------------------------------------------
c     # interpolate_face_ghost
c     # interpolate_corner_ghost
c     # interpolate_to_fine_patch
c     #
c     # Other routines :
c     # compute_slopes (for limited function reconstruction)
c     # fixcapaq (to preserve conservation)
c     #
c     # Note that fixcapaq is only used when regridding;  ghost
c     # cell interpolation is not conservative in the mapped case.
c     # (Should it be?  We are going to correct the flux mixmatch
c     # anyhow, so maybe the accuracy of the ghost cell values is
c     # more important.)
c     # ----------------------------------------------------------


c     # ----------------------------------------------------------
c     # This routine is used for both mapped and non-mapped
c     # cases.
c     # ----------------------------------------------------------
      subroutine fc2d_geoclaw_fort_interpolate_face(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mbathy,
     &      idir,iface_coarse,num_neighbors,refratio,igrid,
     &      transform_ptr)

      use geoclaw_module, ONLY:dry_tolerance, sea_level
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer maux, mbathy
      integer num_neighbors
      integer*8 transform_ptr
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq,r2, m
      integer i, ic1, ic2, ibc, ifine,i1
      integer j, jc1, jc2, jbc, jfine,j1
      integer ic_add, jc_add, ic, jc, mth
      double precision gradx, grady, qc, sl, sr, value
      double precision fclaw2d_clawpatch_compute_slopes
      double precision etabarc(-1:1, -1:1), h, b
      integer iii,jjj

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      logical fclaw2d_clawpatch_is_valid_interp
      logical skip_this_grid

      integer a(2,2), f(2)
      integer ii,jj,dc(2),df(2,0:rr2-1),iff,jff
      double precision shiftx(0:rr2-1),shifty(0:rr2-1)

      mth = 5
      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


      call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)


c     # This needs to be written for refratios .ne. 2.
      m = 0
      do jj = 0,1
         do ii = 0,1
c           # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj
c           # Direction on fine grid (converted using metric). Divide
c           # by refratio to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/refratio
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/refratio

c           # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5d0)/2.d0
            shifty(m) = (jj-0.5d0)/2.d0
            m = m + 1
         enddo
      enddo


c     # Create map :
      do mq = 1,meqn
         if (idir .eq. 0) then

c           # this ensures that we get 'hanging' corners
            if (iface_coarse .eq. 0) then
               ic = 1
            elseif (iface_coarse .eq. 1) then
               ic = mx
            endif

            do jc = 1,mx
               i1 = ic
               j1 = jc
               call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2,
     &               transform_ptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. 
     &        fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  if (mq .eq. 1) then
!                    ! Calculate surface elevation eta using dry limiting
                     do ii = -1, 1
                        do jj = -1, 1
                           h = qcoarse(1,ic+ii,jc+jj)
                           b = aux_coarse(mbathy,ic+ii,jc+jj)
                           if (h < dry_tolerance) then
                              etabarc(ii,jj) = sea_level
                           else
                              etabarc(ii,jj) = h + b
                           endif
                        enddo
                     enddo

                     qc = etabarc(0,0)

                     sl = qc - etabarc(-1,0)
                     sr = etabarc(1,0) - qc
                     gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     sl = qc - etabarc(0,-1)
                     sr = etabarc(0,1) - qc
                     grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)
                     do m = 0,rr2-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m)
     &                        -aux_fine(mbathy,iff,jff)
                        qfine(mq,iff,jff) = max(value,0.0)
                     enddo
                  else
                     qc = qcoarse(mq,ic,jc)

                     sl = (qc - qcoarse(mq,ic-1,jc))
                     sr = (qcoarse(mq,ic+1,jc) - qc)
                     gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     sl = (qc - qcoarse(mq,ic,jc-1))
                     sr = (qcoarse(mq,ic,jc+1) - qc)
                     grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     do m = 0,rr2-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m)
                        qfine(mq,iff,jff) = value
                     enddo
                  endif
c-----------------------------------------------------------------------------
               endif
            enddo
         else                           ! other direction
            if (iface_coarse .eq. 2) then
               jc = 1
            elseif (iface_coarse .eq. 3) then
               jc = my
            endif

            do ic = 1,mx
    1          i1 = ic
               j1 = jc
               call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2,
     &               transform_ptr)
c              # ---------------------------------------------
c              # Two 'half-size' neighbors will be passed into
c              # this routine.  Only half of the coarse grid ghost
c              # indices will be valid for the particular grid
c              # passed in.  We skip those ghost cells that will
c              # have to be filled in by the other half-size
c              # grid.
c              # ---------------------------------------------
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. 
     &        fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo

               if (.not. skip_this_grid) then
!                 Calculate surface elevation eta using dry limiting
                  if (mq .eq. 1) then
                     do ii = -1, 1
                        do jj = -1, 1
                           h = qcoarse(1,ic+ii,jc+jj)
                           b = aux_coarse(mbathy,ic+ii,jc+jj)
                           if (h < dry_tolerance) then
                              etabarc(ii,jj) = sea_level
                           else
                              etabarc(ii,jj) = h + b
                           endif
                        enddo
                     enddo

                     qc = etabarc(0,0)

                     sl = qc - etabarc(-1,0)
                     sr = etabarc(1,0) - qc
                     gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     sl = qc - etabarc(0,-1)
                     sr = etabarc(0,1) - qc
                     grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     do m = 0,rr2-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m)
     &                        - aux_fine(mbathy,iff,jff)
                        qfine(mq,iff,jff) = max(value,0.0)
                     enddo
                  else
                     qc = qcoarse(mq,ic,jc)

                     sl = (qc - qcoarse(mq,ic-1,jc))
                     sr = (qcoarse(mq,ic+1,jc) - qc)
                     gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                     sl = (qc - qcoarse(mq,ic,jc-1))
                     sr = (qcoarse(mq,ic,jc+1) - qc)
                     grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)
                     do m = 0,rr2-1
                        iff = i2(0) + df(1,m)
                        jff = j2(0) + df(2,m)
                        value = qc + gradx*shiftx(m) + grady*shifty(m)
                        qfine(mq,iff,jff) = value
                     enddo
                  endif
               endif                    !! Don't skip this grid
            enddo                       !! i loop
         endif                          !! end idir branch
      enddo                             !! endo mq loop

      end

      subroutine fc2d_geoclaw_fort_interpolate_corner(mx,my,mbc,meqn,
     &      refratio,qcoarse,qfine,maux,aux_coarse,aux_fine,mbathy,
     &      icorner_coarse,transform_ptr)

      use geoclaw_module, ONLY:dry_tolerance, sea_level
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer mbathy,maux
      integer*8 transform_ptr
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ic, jc, mq, ibc,jbc, mth,i,j
      double precision qc, sl, sr, gradx, grady
      double precision fclaw2d_clawpatch_compute_slopes, value

      double precision etabarc(-1:1, -1:1), h, b

c     # This should be refratio*refratio.
      integer i1,j1,m, r2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      integer a(2,2), f(2)
      integer ii,jj,iff,jff,dc(2),df(2,0:rr2-1)
      double precision shiftx(0:rr2-1), shifty(0:rr2-1)
      logical check_indices

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

      call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)

      m = 0
      do jj = 0,1
         do ii = 0,1
c           # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

c           # Direction on fine grid (converted using metric). Divide
c           # by 2 (refratio) to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/2
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/2

c           # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5d0)/2.d0
            shifty(m) = (jj-0.5d0)/2.d0
            m = m + 1
         enddo
      enddo


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

c     # Interpolate coarse grid corners to fine grid corner ghost cells
      i1 = ic
      j1 = jc
      call fclaw2d_clawpatch_transform_corner_half(i1,j1,i2,j2,
     &      transform_ptr)

      do mq = 1,meqn
         if (mq .eq. 1) then
            do ii = -1, 1
               do jj = -1, 1
                  h = qcoarse(1,ic+ii,jc+jj)
                  b = aux_coarse(mbathy,ic+ii,jc+jj)
                  if (h < dry_tolerance) then
                     etabarc(ii,jj) = sea_level
                  else
                     etabarc(ii,jj) = h + b
                  endif
               enddo
            enddo

            qc = etabarc(0,0)

            sl = qc - etabarc(-1,0)
            sr = etabarc(1,0) - qc
            gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

            sl = qc - etabarc(0,-1)
            sr = etabarc(0,1) - qc
            grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

!!            qc = qcoarse(mq,ic,jc) + aux_coarse(mbathy,ic,jc)

            do m = 0,rr2-1
               iff = i2(0) + df(1,m)
               jff = j2(0) + df(2,m)
               value = qc + gradx*shiftx(m) + grady*shifty(m)
     &               - aux_fine(mbathy,iff,jff)
               qfine(mq,iff,jff) = max(value,0.0)
            enddo
         else
            qc = qcoarse(mq,ic,jc)

            sl = (qc - qcoarse(mq,ic-1,jc))
            sr = (qcoarse(mq,ic+1,jc) - qc)
            gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

            sl = (qc - qcoarse(mq,ic,jc-1))
            sr = (qcoarse(mq,ic,jc+1) - qc)
            grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)
            do m = 0,rr2-1
               iff = i2(0) + df(1,m)
               jff = j2(0) + df(2,m)
               value = qc + gradx*shiftx(m) + grady*shifty(m)
               qfine(mq,iff,jff) = value
            enddo
         endif

      enddo

      end

      subroutine fc2d_geoclaw_fort_interpolate2fine(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mbathy,
     &      p4est_refineFactor,refratio,igrid)

      use geoclaw_module, ONLY:dry_tolerance, sea_level
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid, maux,mbathy
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      integer ic,jc,ic_add, jc_add, iff, jf
      double precision qc, qf, shiftx, shifty, sl, sr, gradx, grady
      double precision fclaw2d_clawpatch_compute_slopes
      double precision uc(-1:1,-1:1), uf
      double precision etabarc(-1:1, -1:1), h, b, u, hfsum
      double precision coarseumin, coarseumax
      logical redefine

      hfsum = 0.d0
c     # Use minmod to calculate slope.
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
c        # First loop over quadrant (i1,i2)x(j1,j2) of the coarse grid
         do i = i1,i2
            do j = j1,j2
               ic = i + ic_add
               jc = j + jc_add
               ! Calculate surface elevation eta using dry limiting
               do ii = -1, 1
                  do jj = -1, 1
                     h = qcoarse(1,ic+ii,jc+jj)
                     b = aux_coarse(mbathy,ic+ii,jc+jj)
                     if (h < dry_tolerance) then
                        etabarc(ii,jj) = sea_level
                     else
                        etabarc(ii,jj) = h + b
                     endif
                  enddo
               enddo
c--------------start interpolation
               if (mq .eq. 1) then
c                 # Interpolate sea surface height rather than just the
c                 # water column height.
                  qc = etabarc(0,0)

                  sl = qc - etabarc(-1,0)
                  sr = etabarc(1,0) - qc
                  gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                  sl = qc - etabarc(0,-1)
                  sr = etabarc(0,1) - qc
                  grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

c                 # Fill in fine grid values from coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qf = qc + shiftx*gradx + shifty*grady -
     &                        aux_fine(mbathy,iff,jf)
                        qfine(mq,iff,jf) = qf
                     enddo
                  enddo
               else
c                 # interpolate momentum components in the usual way.
c                 # But then make sure that no new extrema are created.
                  qc = qcoarse(mq,ic,jc)

                  sl = (qc - qcoarse(mq,ic-1,jc))
                  sr = (qcoarse(mq,ic+1,jc) - qc)
                  gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(mq,ic,jc-1))
                  sr = (qcoarse(mq,ic,jc+1) - qc)
                  grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

c                 # Fill in refined values on coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qfine(mq,iff,jf) =
     &                        qc + shiftx*gradx + shifty*grady
                        hfsum = hfsum + qfine(mq,iff,jf)
                     enddo
                  enddo


c------------check to make sure we are not creating any new extrema

c                 # calculate coarse cells' velocity
c                  do ii = -1,1
c                     do jj = -1,1
c                        if (qcoarse(1,ic+ii,jc+jj) .eq. 0.d0) then
c                           uc(ii,jj) = 0.d0
c                        else
c                           uc(ii,jj) = qcoarse(mq,ic+ii,jc+jj)/
c     &                                 qcoarse(1,ic+ii,jc+jj)
c                        endif
c                     enddo
c                  enddo
c                 # find the maximum/minimum velocities among coarse cells
c                  coarseumax = -1d99
c                  coarseumin = 1d99
c                  do ii = -1,1
c                     do jj = -1,1
c                        coarseumax = max(coarseumax,uc(ii,jj))
c                        coarseumin = min(coarseumin,uc(ii,jj))
c                     enddo
c                  enddo

c                  redefine = .false.
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        iff = (i-1)*refratio + ii
c                        jf = (j-1)*refratio + jj
c                        if (qfine(1,iff,jf) .eq. 0.d0) then
c                           uf = 0.d0
c                        else
c                           uf = qfine(mq,iff,jf)/qfine(1,iff,jf)
c                        endif
c                        if (uf .gt. coarseumax .or. uf .lt. coarseumin)
c     &                        then
c                           redefine = .true.
c                        endif
c                     enddo
c                  enddo

c                  if (redefine) then
c                     u = qcoarse(mq,ic,jc)/
c     &                 max(qcoarse(1,ic,jc), hfsum/(refratio*refratio))
c                     do ii = 1,refratio
c                        do jj = 1,refratio
c                           iff = (i-1)*refratio + ii
c                           jf = (j-1)*refratio + jj
c                           qfine(mq,iff,jf) = qfine(1,iff,jf)*u
c                        enddo
c                     enddo
c                  endif
c------end of checking to make sure we are not creating any new extrema

               endif
c------end of interpolation
            enddo
         enddo
      enddo

      end