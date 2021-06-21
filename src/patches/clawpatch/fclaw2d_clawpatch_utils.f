      double precision function 
     &      fclaw2d_clawpatch_compute_slopes(sl,sr,mth)
      implicit none

      double precision sl,sr, sc, philim, slim
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
         fclaw2d_clawpatch_compute_slopes = slim*sl
      else
c        # Use AMRClaw slopes  (use minimum in absolute value;  sign is
c        # chosen from centered (sc) slope
         sc = (sl + sr)/2.d0
         fclaw2d_clawpatch_compute_slopes = 
     &         min(2*abs(sl),2*abs(sr),abs(sc))*
     &         max(0.d0,sign(1.d0,sl*sr))*sign(1.d0,sc)

c        # Do this to guarantee that ghost cells are used; this is a check
c        # on the ghost-fill procedures.  Could raise an exception if face
c        # a patch communicates with more two or more procs.  If this
c        # is uncommented, also uncomment warning in fclaw2d_ghost_fill.cpp
c         fclaw2d_clawpatch_compute_slopes = sc

      endif


      end

      subroutine fclaw2d_clawpatch_build_transform(transform_ptr,a,f)
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
      call fclaw2d_clawpatch_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      f(1) = mi(1)
      f(2) = mj(1)

      i1 = 1
      j1 = 0
      call fclaw2d_clawpatch_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      a(1,1) = mi(1) - f(1)
      a(2,1) = mj(1) - f(2)

      i1 = 0
      j1 = 1
      call fclaw2d_clawpatch_transform_face_half(i1,j1,mi,mj,
     &      transform_ptr)
      a(1,2) = mi(1) - f(1)
      a(2,2) = mj(1) - f(2)

      end

      logical function fclaw2d_clawpatch_check_indices(iff,jff,i2,j2)
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

      fclaw2d_clawpatch_check_indices = found_iff .and. found_jff


      end

      logical function fclaw2d_clawpatch_is_valid_interp(i,j,mx,my,mbc)
      implicit none
      integer i,j,mx, my, mbc

      logical i1, j1

      i1 = 1-mbc .le. i .and. i .le. mx+mbc
      j1 = 1-mbc .le. j .and. j .le. my+mbc

      fclaw2d_clawpatch_is_valid_interp = i1 .and. j1

      end

      logical function fclaw2d_clawpatch_is_valid_average(i,j,mx,my)
      implicit none

      integer i,j,mx,my
      logical i1, j1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my

      fclaw2d_clawpatch_is_valid_average = i1 .and. j1

      end



