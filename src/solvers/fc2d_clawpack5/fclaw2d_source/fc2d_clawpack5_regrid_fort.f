c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

      subroutine fc2d_clawpack5_fort_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision qmin, qmax

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1)
      qmax = q(mq,1,1)
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            if (qmax - qmin .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end


c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine fc2d_clawpack5_fort_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(mq,1,1)
      qmax = q0(mq,1,1)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call fc2d_clawpack5_get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call fc2d_clawpack5_get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call fc2d_clawpack5_get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call fc2d_clawpack5_get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      coarsen_threshold,tag_patch)

      end

      subroutine fc2d_clawpack5_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer i,j

      do i = 1,mx
         do j = 1,my
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            if (qmax - qmin .gt. coarsen_threshold) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end


c     # Conservative intepolation to fine grid patch
      subroutine fc2d_clawpack5_fort_interpolate2fine(mx,my,mbc,meqn,
     &      qcoarse, qfine, areacoarse, areafine, igrid, manifold)
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
      double precision compute_slopes

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

      if (manifold .ne. 0) then
         call fc2d_clawpack5_fort_fixcapaq2(mx,my,mbc,meqn,
     &         qcoarse,qfine,areacoarse,areafine,igrid)
      endif


      end

c> \ingroup  Averaging
c> Average fine grid siblings to parent coarse grid.
      subroutine fc2d_clawpack5_fort_average2coarse(mx,my,mbc,meqn,
     &      qcoarse,qfine,areacoarse,areafine,igrid,manifold)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      integer manifold
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      integer mq
      double precision sum
      logical is_manifold

c     # This should be refratio*refratio.
      integer i1,j1, r2, m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      double precision kc, kf, qf

      p4est_refineFactor = 2
      refratio = 2

      is_manifold = manifold .eq. 1

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


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
               i1 = i+ic_add
               j1 = j+jc_add
               m = 0
               do jj = 1,refratio
                  do ii = 1,refratio
                     i2(m) = (i-1)*refratio + ii
                     j2(m) = (j-1)*refratio + jj
                     m = m + 1
                  enddo
               enddo
               if (is_manifold) then
                  sum = 0
                  do m = 0,r2-1
                     qf = qfine(mq,i2(m),j2(m))
                     kf = areafine(i2(m),j2(m))
                     sum = sum + kf*qf
                  enddo
                  kc = areacoarse(i1,j1)
                  qcoarse(mq,i1,j1) = sum/kc
               else
                  sum = 0
                  do m = 0,r2-1
                     qf = qfine(mq,i2(m),j2(m))
                     sum = sum + qf
                  enddo
                  qcoarse(mq,i1,j1) = sum/r2
               endif
            enddo
         enddo
      enddo
      end

c     # ------------------------------------------------------
c     # So far, this is only used by the interpolation from
c     # coarse to fine when regridding.  But maybe it should
c     # be used by the ghost cell routines as well?
c     # ------------------------------------------------------
      subroutine fc2d_clawpack5_fort_fixcapaq2(mx,my,mbc,meqn,
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
