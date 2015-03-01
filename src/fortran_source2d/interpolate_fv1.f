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
     &      transform_ptr)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      integer*8 transform_ptr
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2, m
      integer i, ic1, ic2, ibc, ii, ifine
      integer j, jc1, jc2, jbc, jj, jfine
      integer ic_add, jc_add, ic, jc, mth
      double precision shiftx, shifty, gradx, grady, qc, sl, sr, value
      double precision compute_slopes

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      logical is_valid_interp
      logical skip_this_grid

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

            if (iface_coarse .eq. 0) then
               ic = 1
            elseif (iface_coarse .eq. 1) then
               ic = mx
            endif
            do jc = 1,mx
               call fclaw2d_transform_face_half(ic,jc,i2,j2,
     &               transform_ptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  qc = qcoarse(ic,jc,mq)
c                 # Compute limited slopes in both x and y. Note we are not
c                 # really computing slopes, but rather just differences.
c                 # Scaling is accounted for in 'shiftx' and 'shifty', below.
                  sl = (qc - qcoarse(ic-1,jc,mq))
                  sr = (qcoarse(ic+1,jc,mq) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(ic,jc-1,mq))
                  sr = (qcoarse(ic,jc+1,mq) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # This works for smooth grid mappings as well.
                  do m = 0,r2-1
                     shiftx =
     &                     (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                     shifty =
     &                     (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                     value = qc + shiftx*gradx + shifty*grady
                     qfine(i2(m),j2(m),mq) = value
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
               call fclaw2d_transform_face_half(ic,jc,i2,j2,
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
                  if (.not. is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  qc = qcoarse(ic,jc,mq)

                  sl = (qc - qcoarse(ic-1,jc,mq))
                  sr = (qcoarse(ic+1,jc,mq) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(ic,jc-1,mq))
                  sr = (qcoarse(ic,jc+1,mq) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # Compute interpolant for each of four
c                 # fine grid cells.
                  do m = 0,r2-1
                     shiftx =
     &                     (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
                     shifty =
     &                     (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
                     value = qc + shiftx*gradx + shifty*grady
                     qfine(i2(m),j2(m),mq) = value
                  enddo
               endif
            enddo
         endif
      enddo

      end

      logical function is_valid_interp(i,j,mx,my,mbc)
      implicit none
      integer i,j,mx, my, mbc

      logical i1, i2
      logical j1, j2

      i1 = 1-mbc .le. i .and. i .le. mx+mbc
      j1 = 1-mbc .le. j .and. j .le. my+mbc

      is_valid_interp = i1 .and. j1

      end

      subroutine interpolate_corner_ghost(mx,my,mbc,meqn,
     &      refratio,
     &      qcoarse,qfine,icorner_coarse,transform_ptr)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer*8 transform_ptr
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
         call fclaw2d_transform_corner_half(i1,j1,i2,j2,
     &         transform_ptr)
         do m = 0,r2-1
            shiftx = (i2(m)-i2(0)- refratio/2.d0 + 0.5)/refratio
            shifty = (j2(m)-j2(0)- refratio/2.d0 + 0.5)/refratio
            value = qc + shiftx*gradx + shifty*grady
            qfine(i2(m),j2(m),mq) = value
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
c        # Use AMRClaw slopes  (use minimum in absolute value;  sign is
c        # chosen from centered (sc) slope
         sc = (sl + sr)/2.d0
         compute_slopes = min(abs(sl),abs(sr),abs(sc))*
     &         max(0.d0,sign(1.d0,sl*sr))*sign(1.d0,sc)
      endif


      end
