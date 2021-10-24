c> @file
c> Clawpatch utilities

c -------------------------------------------------------
c> @brief Slope limiting function
c>
c> @param[in] sl, sr the slope on the left and right sides
c> @param[in] mth type of limiting to use
c> @return the limited slope
c -------------------------------------------------------
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

c     ---------------------------------------------------------------
c>    @brief Determines the tranformation from a coarse patch to it's
c>    half-sized neighbor's coordinate system
c>
c>    @param[in]  the pointer to the fclaw2d_patch_transform_data struct
c>    @param[out] a the 2x2 tranform matrix for the i,j indexes of a patch
c>    @param[out] f the transform vector
c     ---------------------------------------------------------------
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

c --------------------------------------------------------------------
c> @brief checks if the index is in the range given
c>
c> @param[in] iff, jff the index to check
c> @param[in] i2, j2 the list of indecies
c> @return true if the index is valid
c --------------------------------------------------------------------
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

c --------------------------------------------------------------------
c> @brief checks if the index is valid for interpolation
c>
c> @param[in] i, j the idnex to check
c> @param[in] mx, my the number of cells in the x and y directions
c> @param[in] mbc the number of ghost cells
c> @return true if the index is valid
c --------------------------------------------------------------------
      logical function fclaw2d_clawpatch_is_valid_interp(i,j,mx,my,mbc)
      implicit none
      integer i,j,mx, my, mbc

      logical i1, j1

      i1 = 1-mbc .le. i .and. i .le. mx+mbc
      j1 = 1-mbc .le. j .and. j .le. my+mbc

      fclaw2d_clawpatch_is_valid_interp = i1 .and. j1

      end

c --------------------------------------------------------------------
c> @brief checks if the index is valid for averaging
c>
c> @param[in] i, j the idnex to check
c> @param[in] mx, my the number of cells in the x and y directions
c> @return true if the index is valid
c --------------------------------------------------------------------
      logical function fclaw2d_clawpatch_is_valid_average(i,j,mx,my)
      implicit none

      integer i,j,mx,my
      logical i1, j1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my

      fclaw2d_clawpatch_is_valid_average = i1 .and. j1

      end

c --------------------------------------------------------------------
c> Initializes the area and edge length arrays for fclaw2d_clawpatch_registers
c>
c> @param[in]  mx, my number of cells in the x and y directions
c> @param[in]  mbc number of ghost cells
c> @param[in]  dx, dy spacings in the x and y direcitons
c> @param[in]  area area of each cell
c> @param[in]  edgelengths edge lengts of each cell along the x and y directions
c> @param[out] area0, area1, area2, area3 areas of cells along the edges
c> @param[out] el0, el1, el2, el3 edges lengths of the cells along the dges
c> @param[in]  manifold 1 if using manifold
c --------------------------------------------------------------------
      subroutine clawpatch_time_sync_setup(mx,my,mbc,dx,dy,
     &                                     area,edgelengths,
     &                                     area0,area1,area2,area3,
     &                                     el0, el1, el2, el3, 
     &                                     manifold)


      implicit none

      integer mx,my,mbc
      integer manifold
      double precision dx,dy

      double precision  area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision  edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)

      double precision area0(my), area1(my), area2(mx), area3(mx)
      double precision el0(my), el1(my), el2(mx), el3(mx)

      integer i,j
      double precision dxdy

c      include "fclaw2d_metric_terms.i"

      if (manifold .eq. 1) then
         do j = 1,my
            area0(j) = area(1,j)
            area1(j) = area(mx,j)
            el0(j) = edgelengths(1,j,1)
            el1(j) = edgelengths(mx+1,j,1)
         enddo
    
         do i = 1,mx
            area2(i) = area(i,1)
            area3(i) = area(i,my)
            el2(i) = edgelengths(i,1,2)
            el3(i) = edgelengths(i,my+1,2)
         enddo
      else
         dxdy = dx*dy
         do j = 1,my
            area0(j) = dxdy
            area1(j) = dxdy
            el0(j) = dy
            el1(j) = dy
         enddo
         do i = 1,mx
            area2(i) = dxdy
            area3(i) = dxdy
            el2(i) = dx
            el3(i) = dx
         enddo
      endif


      end



