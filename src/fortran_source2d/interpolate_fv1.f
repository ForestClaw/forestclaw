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
      subroutine interpolate_face_ghost(mx,my,mbc,meqn,
     &      qcoarse,qfine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid,
     &      transform_cptr)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      integer*8 transform_cptr
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2, m
      integer i, ic1, ic2, ibc, ii, ifine
      integer j, jc1, jc2, jbc, jj, jfine
      integer ic_add, jc_add, ic, jc, mth
      double precision shiftx, shifty, gradx, grady, qc, sl, sr, value
      double precision compute_slopes

c     # This should be refratio*refratio.
      integer i1,j1
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      mth = 5
      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

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

c              # This works for smooth grid mappings as well.
               i1 = ic
               j1 = jc
               call transform_face_halfsize(i1,j1,i2,j2,transform_cptr)
               do m = 0,r2-1
                  shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                  shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                  value = qc + shiftx*gradx + shifty*grady
                  qfine(i2(m),j2(m),mq) = value
               enddo
            enddo
         else
            ic_add = igrid*mx/num_neighbors
c           # this ensures that we get 'hanging' corners
            ic1 = 1-igrid
            ic2 = mx/num_neighbors + (1-igrid)

            if (iface_coarse .eq. 2) then
               jc = 1
            elseif (iface_coarse .eq. 3) then
               jc = my
            endif
            do i = ic1, ic2
               ic = i + ic_add
               qc = qcoarse(ic,jc,mq)

               sl = (qc - qcoarse(ic-1,jc,mq))
               sr = (qcoarse(ic+1,jc,mq) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(ic,jc-1,mq))
               sr = (qcoarse(ic,jc+1,mq) - qc)
               grady = compute_slopes(sl,sr,mth)

c              # Compute interpolant for each of four
c              # fine grid cells.
               i1 = ic
               j1 = jc
               call transform_face_halfsize(i1,j1,i2,j2,transform_cptr)
               do m = 0,r2-1
                  shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                  shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                  value = qc + shiftx*gradx + shifty*grady
                  qfine(i2(m),j2(m),mq) = value
               enddo
            enddo
         endif
      enddo

      end

      subroutine interpolate_corner_ghost(mx,my,mbc,meqn,
     &      refratio,
     &      qcoarse,qfine,icorner_coarse,transform_cptr)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer*8 transform_cptr
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ic, jc, mq, ibc,jbc, mth,i,j
      double precision qc, sl, sr, gradx, grady, shiftx, shifty
      double precision compute_slopes, value

c     # This should be refratio*refratio.
      integer i1,j1,m, r2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


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

c        # This works for smooth grid mappings as well.
         i1 = ic
         j1 = jc
         call transform_corner_halfsize2(i1,j1,i2,j2,transform_cptr)
         do m = 0,r2-1
            shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
            shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
            value = qc + shiftx*gradx + shifty*grady
            qfine(i2(m),j2(m),mq) = value
         enddo
      enddo

      end


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

c     # ------------------------------------------------------
c     # So far, this is only used by the interpolation from
c     # coarse to fine when regridding.  But maybe it should
c     # be used by the ghost cell routines as well?
c     # ------------------------------------------------------
      subroutine fixcapaq2(mx,my,mbc,meqn,
     &      qcoarse,qfine,
     &      areacoarse,areafine,
     &      p4est_refineFactor,
     &      refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn, refratio, igrid
      integer p4est_refineFactor

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
      double precision kf, kc, r2, sum, cons_diff, qf, qc

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
                     qf = qfine(ifine,jfine,m)
                     sum = sum + kf*qf
                  enddo
               enddo

               kc = areacoarse(i+ic_add,j+jc_add)
               qc = qcoarse(i+ic_add, j+jc_add,m)
               cons_diff = (qc*kc - sum)/r2

               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine  = (i-1)*refratio + ii
                     jfine  = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     if (kf .eq. 0) then
                        write(6,*) 'fixcapaq : kf = 0'
                        stop
                     endif
                     qfine(ifine,jfine,m) = qfine(ifine,jfine,m) +
     &                     cons_diff/kf
                  enddo
               enddo
            enddo  !! end of meqn
         enddo
      enddo

      end
