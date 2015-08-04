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
      integer i, ic1, ic2, ibc, ifine,i1
      integer j, jc1, jc2, jbc, jfine,j1
      integer ic_add, jc_add, ic, jc, mth
      double precision gradx, grady, qc, sl, sr, value
      double precision compute_slopes

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      logical is_valid_interp
      logical skip_this_grid

      integer a(2,2)
      integer ii,jj,dc(2),df(2,0:rr2-1),iff,jff
      double precision shiftx(0:rr2-1),shifty(0:rr2-1)

      mth = 5
      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


      call build_transform(transform_ptr,a)

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
               call fclaw2d_transform_face_half(i1,j1,i2,j2,
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

                  do m = 0,rr2-1
                     iff = i2(0) + df(1,m)
                     jff = j2(0) + df(2,m)
                     value = qc + gradx*shiftx(m) + grady*shifty(m)
                     qfine(iff,jff,mq) = value
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
    1          i1 = ic
               j1 = jc
               call fclaw2d_transform_face_half(i1,j1,i2,j2,
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

                  do m = 0,rr2-1
                     iff = i2(0) + df(1,m)
                     jff = j2(0) + df(2,m)
                     value = qc + gradx*shiftx(m) + grady*shifty(m)
                     qfine(iff,jff,mq) = value
                  enddo

               endif                    !! Don't skip this grid
            enddo                       !! i loop
         endif                          !! end idir branch
      enddo                             !! endo mq loop

      end

      logical function check_indices(iff,jff,i2,j2)
      implicit none

      integer iff,jff,i2(0:3),j2(0:3)
      integer m
      logical found_iff, found_jff

      found_iff = .false.
      do m = 0,3
         if (i2(m) .eq. iff) then
            found_iff = .true.
            exit
         endif
      enddo

      found_jff = .false.
      do m = 0,3
         if (j2(m) .eq. jff) then
            found_jff = .true.
            exit
         endif
      enddo

      check_indices = found_iff .and. found_jff


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
      double precision qc, sl, sr, gradx, grady
      double precision compute_slopes, value

c     # This should be refratio*refratio.
      integer i1,j1,m, r2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      integer a(2,2)
      integer ii,jj,iff,jff,dc(2),df(2,0:rr2-1)
      double precision shiftx(0:rr2-1), shifty(0:rr2-1)
      logical check_indices

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

      call build_transform(transform_ptr,a)

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

      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)

c        # Interpolate coarse grid corners to fine grid corner ghost cells
         i1 = ic
         j1 = jc
         call fclaw2d_transform_face_half(i1,j1,i2,j2,
     &         transform_ptr)


c        # Compute limited slopes in both x and y. Note we are not
c        # really computing slopes, but rather just differences.
c        # Scaling is accounted for in 'shiftx' and 'shifty', below.
         sl = (qc - qcoarse(ic-1,jc,mq))
         sr = (qcoarse(ic+1,jc,mq) - qc)
         gradx = compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(ic,jc-1,mq))
         sr = (qcoarse(ic,jc+1,mq) - qc)
         grady = compute_slopes(sl,sr,mth)

         do m = 0,rr2-1
            iff = i2(0) + df(1,m)
            jff = j2(0) + df(2,m)
            value = qc + gradx*shiftx(m) + grady*shifty(m)
            qfine(iff,jff,mq) = value
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

c        # Do this to guarantee that ghost cells are used; this is a check
c        # on the ghost-fill procedures.  Could raise an exception if face
c        # a patch communicates with more two or more procs.  If this
c        # is uncommented, also uncomment warning in fclaw2d_ghost_fill.cpp
c         compute_slopes = sc

      endif


      end

      subroutine build_transform(transform_ptr,a)
      implicit none

      integer a(2,2)
      integer*8 transform_ptr
      integer f(2)
      integer mi(4),mj(4)
      integer i1,j1

c     # Assume index mapping fclaw2d_transform_face_half has the
c     # the form
c     #
c     #       T(ic,jc) = A*(ic,jc) + F = (iff,jff)
c     #
c     # where (ic,jc) is the coarse grid index, and (iff,jff)
c     # is the first fine grid index.
c     #
c     # We can recover coefficients A(2,2) with the following
c     # calls to T.

      i1 = 0
      j1 = 0
      call fclaw2d_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      f(1) = mi(1)
      f(2) = mj(1)

      i1 = 1
      j1 = 0
      call fclaw2d_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      a(1,1) = mi(1) - f(1)
      a(2,1) = mj(1) - f(2)

      i1 = 0
      j1 = 1
      call fclaw2d_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      a(1,2) = mi(1) - f(1)
      a(2,2) = mj(1) - f(2)

      end
