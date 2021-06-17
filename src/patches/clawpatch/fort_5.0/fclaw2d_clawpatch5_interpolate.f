c     # ----------------------------------------------------------
c     # Interpolation routines - (i,j,mq) ordering
c     # ----------------------------------------------------------
c     # interpolate_face_ghost
c     # interpolate_corner_ghost
c     # interpolate_to_fine_patch
c     #
c     # Other routines :
c     # fclaw2d_clawpatch_compute_slopes (for limited function reconstruction)
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
      subroutine fclaw2d_clawpatch5_fort_interpolate_face(mx,my,mbc,
     &      meqn, qcoarse,qfine, idir,iface_coarse,num_neighbors,
     & refratio,igrid,transform_ptr)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      integer*8 transform_ptr
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq,r2, m
      integer i1
      integer j1
      integer ic, jc, mth
      double precision gradx, grady, qc, sl, sr, value
      double precision fclaw2d_clawpatch_compute_slopes

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
     &       fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  qc = qcoarse(mq,ic,jc)
c                 # Compute limited slopes in both x and y. Note we are not
c                 # really computing slopes, but rather just differences.
c                 # Scaling is accounted for in 'shiftx' and 'shifty', below.
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
         else
            if (iface_coarse .eq. 2) then
               jc = 1
            elseif (iface_coarse .eq. 3) then
               jc = my
            endif
            do ic = 1,mx
               i1 = ic
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
     &       fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
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

               endif                    !! Don't skip this grid
            enddo                       !! i loop
         endif                          !! end idir branch
      enddo                             !! endo mq loop

      end

      subroutine fclaw2d_clawpatch5_fort_interpolate_corner(mx,my,
     &      mbc,meqn,refratio, qcoarse,qfine,icorner_coarse,
     &      transform_ptr)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer*8 transform_ptr
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ic, jc, mq, mth
      double precision qc, sl, sr, gradx, grady
      double precision fclaw2d_clawpatch_compute_slopes, value

c     # This should be refratio*refratio.
      integer i1,j1,m, r2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      integer a(2,2), f(2)
      integer ii,jj,iff,jff,dc(2),df(2,0:rr2-1)
      double precision shiftx(0:rr2-1), shifty(0:rr2-1)

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
         qc = qcoarse(mq,ic,jc)

c        # Compute limited slopes in both x and y. Note we are not
c        # really computing slopes, but rather just differences.
c        # Scaling is accounted for in 'shiftx' and 'shifty', below.
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

      enddo

      end

c     # Conservative intepolation to fine grid patch
      subroutine fclaw2d_clawpatch5_fort_interpolate2fine(mx,my,mbc,
     &    meqn, qcoarse, qfine, areacoarse, areafine, igrid, manifold)
      implicit none

      integer mx,my,mbc,meqn
      integer igrid, manifold

      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      integer ic,jc,ic_add, jc_add
      double precision qc, shiftx, shifty, sl, sr, gradx, grady
      double precision fclaw2d_clawpatch_compute_slopes

      integer p4est_refineFactor,refratio

      p4est_refineFactor = 2
      refratio = 2

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
               gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(mq,ic,jc-1))
               sr = (qcoarse(mq,ic,jc+1) - qc)
               grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

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

      if (manifold .ne. 0) then
         call fclaw2d_clawpatch5_fort_fixcapaq2(mx,my,mbc,meqn,
     &         qcoarse,qfine,areacoarse,areafine,igrid)
      endif


      end


c     # ------------------------------------------------------
c     # So far, this is only used by the interpolation from
c     # coarse to fine when regridding.  But maybe it should
c     # be used by the ghost cell routines as well?
c     # ------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_fixcapaq2(mx,my,mbc,meqn,
     &      qcoarse,qfine, areacoarse,areafine,igrid)
      implicit none

      integer mx,my,mbc,meqn, refratio, igrid
      integer p4est_refineFactor

      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
      double precision kf, kc, r2, sum, cons_diff, qf, qc

      p4est_refineFactor = 2
      refratio = 2

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*my/p4est_refineFactor

c     # ------------------------------------------------------
c     # This routine ensures that the interpolated solution
c     # has the same mass as the coarse grid solution
c     # -------------------------------------------------------


      r2 = refratio*refratio
      do m = 1,meqn
         do i = 1,mx/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
               sum = 0.d0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     qf = qfine(m,ifine,jfine)
                     sum = sum + kf*qf
                  enddo
               enddo

               kc = areacoarse(i+ic_add,j+jc_add)
               qc = qcoarse(m,i+ic_add, j+jc_add)
               cons_diff = (qc*kc - sum)/r2

               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine  = (i-1)*refratio + ii
                     jfine  = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     qfine(m,ifine,jfine) = qfine(m,ifine,jfine) +
     &                     cons_diff/kf
                  enddo
               enddo
            enddo  !! end of meqn
         enddo
      enddo

      end
