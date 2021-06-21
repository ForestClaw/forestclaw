c> \file
c> \defgroup Averaging Average fine grids to a coarse grid
c> Average cells from coarse grid to fine grid.
c>
c> Routines described here average are used to fill coarse grid ghost
c> cells, and average sibling grids onto a parent grid.  Indices
c> for cells at block boundaries are transformed using encodings
c> stored in `transform_cptr`.
c>
c> \param [in] mx,my       Number of cells in x,y direction
c> \param [in] mbc      Number of ghost cells
c> \param [in] meqn     Number of equations
c> \param [in] qcoarse,qfine  Solution on coarse,fine grid
c> \param [in] areacoarse,areafine  Area of mesh cells on coarse,fine grids.
c> \param [in] idir     Face orientation - 0 for x-faces; 1 for y-faces [0-1]
c> \param [in] iface    Face number of fine grid [0-3].
c> \param [in] iface_coarse Face number of coarse grid [0-3].
c> \param [in] num_neighbors Number of fine grid neighbors [2].
c> \param [in] refratio  Refinement ratio between coarse and fine grids [2].
c> \param [in] manifold  Flag indicating whether we are on mapped grid [0-1].
c> \param [in] transform_cptr  Encoding for indices at block boundaries (C only).

c> \ingroup Averaging
c> Average fine ghost cell values.
c>
c> Average fine grid interior values to neighboring ghost cell values of
c> the coarse grid.
      subroutine fclaw2d_clawpatch46_fort_average_face(mx,my,mbc,meqn,
     &      qcoarse,qfine,areacoarse, areafine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid,
     &      manifold, transform_cptr)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer manifold
      integer*8 transform_cptr
      integer num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision sum, qf, kf
      logical is_manifold

      integer mq,r2, m
      integer ic, ibc
      integer jc, jbc

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      double precision kc

      logical fclaw2d_clawpatch_is_valid_average, skip_this_grid
      double precision af_sum, qv(0:rr2-1)

      is_manifold = manifold .eq. 1

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

c     # Average fine grid onto coarse grid
      if (idir .eq. 0) then
         do mq = 1,meqn
            do jc = 1,my
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer

                  if (iface_coarse .eq. 0) then
                     ic = 1-ibc
                  elseif (iface_coarse .eq. 1) then
                     ic = mx+ibc
                  endif

                  call 
     &       fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &                  transform_cptr)
c                 # ---------------------------------------------
c                 # Two 'half-size' neighbors will be passed into
c                 # this routine.  Only half of the coarse grid ghost
c                 # indices will be valid for the particular grid
c                 # passed in.  We skip those ghost cells that will
c                 # have to be filled in by the other half-size
c                 # grid.
c                 # ---------------------------------------------
                  skip_this_grid = .false.
                  do m = 0,r2-1
                     if (.not. 
     &       fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my))
     &                     then
                        skip_this_grid = .true.
                        exit
                     endif
                  enddo

                  if (.not. skip_this_grid) then
                     if (is_manifold) then
                        sum = 0
                        af_sum = 0
                        do m = 0,r2-1
                           qf = qfine(i2(m),j2(m),mq)
                           qv(m) = qf
                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                           af_sum = af_sum + kf
                        enddo
c                       # ----------------------------------------
c                       # At block seams, the coarse grid mesh cell
c                       # areas may not have been computed using
c                       # the correct metrics.
c                       # ----------------------------------------
c                       kc = areacoarse(ic,jc)
c                       qcoarse(ic,jc,mq) = sum/kc

c                       # Use areas of the fine grid mesh cells instead.
                        qcoarse(ic,jc,mq) = sum/af_sum
                     else
                        sum = 0
                        do m = 0,r2-1
                           sum = sum + qfine(i2(m),j2(m),mq)
                        enddo
                        qcoarse(ic,jc,mq) = sum/dble(r2)
                     endif
                  endif
               enddo
            enddo
         enddo
      else
         do mq = 1,meqn
c           # idir = 1 (faces 2,3)
            do jbc = 1,mbc
               do ic = 1,mx

                  if (iface_coarse .eq. 2) then
                     jc = 1-jbc
                  elseif (iface_coarse .eq. 3) then
                     jc = my+jbc
                  endif

                  call 
     &      fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &                  transform_cptr)
                  skip_this_grid = .false.
                  do m = 0,r2-1
                     if (.not. 
     &      fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my))
     &                     then
                        skip_this_grid = .true.
                     endif
                  enddo
                  if (.not. skip_this_grid) then
                     if (is_manifold) then
                        sum = 0
                        af_sum = 0
                        do m = 0,r2-1
                           qf = qfine(i2(m),j2(m),mq)
                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                           af_sum = af_sum + kf
                        enddo
                        kc = areacoarse(ic,jc)
c                       qcoarse(ic,jc,mq) = sum/kc
                        qcoarse(ic,jc,mq) = sum/af_sum
                     else
                        sum = 0
                        do m = 0,r2-1
                           sum = sum + qfine(i2(m),j2(m),mq)
                        enddo
                        qcoarse(ic,jc,mq) = sum/dble(r2)
                     endif              !! manifold loop
                  endif                 !! skip grid loop
               enddo
            enddo
         enddo
      endif

      end


c> \ingroup Averaging
c> Average across corners.
      subroutine fclaw2d_clawpatch46_fort_average_corner(mx,my,mbc,meqn,
     &      refratio,qcoarse,qfine,areacoarse,areafine,
     &      manifold,icorner_coarse,transform_cptr)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner_coarse, manifold
      integer*8 transform_cptr
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision sum

      integer ibc,jbc,mq,r2
      logical is_manifold
      double precision qf,kf, kc

c     # This should be refratio*refratio.
      integer i1,j1,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      double precision af_sum

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

      is_manifold = manifold .eq. 1

      r2 = refratio*refratio
c     # Loop over four corner cells on coarse grid
      do ibc = 1,mbc
         do jbc = 1,mbc
c           # Average fine grid corners onto coarse grid ghost corners
            if (icorner_coarse .eq. 0) then
               i1 = 1-ibc
               j1 = 1-jbc
            elseif (icorner_coarse .eq. 1) then
               i1 = mx+ibc
               j1 = 1-jbc
            elseif (icorner_coarse .eq. 2) then
               i1 = 1-ibc
               j1 = my+jbc
            elseif (icorner_coarse .eq. 3) then
               i1 = mx+ibc
               j1 = my+jbc
            endif

c           # Again, a fake routine until the real one is
c           # available (be sure to pass in (i1,j1)
            call fclaw2d_clawpatch_transform_corner_half(i1,j1,i2,j2,
     &            transform_cptr)
            if (is_manifold) then
               do mq = 1,meqn
                  sum = 0
                  af_sum = 0
                  do m = 0,r2-1
                     qf = qfine(i2(m),j2(m),mq)
                     kf = areafine(i2(m),j2(m))
                     sum = sum + kf*qf
                     af_sum = af_sum + kf
                  enddo
                  kc = areacoarse(i1,j1)
c                 qcoarse(i1,j1,mq) = sum/kc
                  qcoarse(i1,j1,mq) = sum/af_sum
               enddo
            else
               do mq = 1,meqn
                  sum = 0
                  do m = 0,r2-1
                     sum = sum + qfine(i2(m),j2(m),mq)
                  enddo
                  qcoarse(i1,j1,mq) = sum/dble(r2)
               enddo
            endif
         enddo
      enddo

      end


c> \ingroup  Averaging
c> Average fine grid siblings to parent coarse grid.
      subroutine fclaw2d_clawpatch46_fort_average2coarse(mx,my,mbc,meqn,
     &      qcoarse,qfine, areacoarse, areafine, igrid,manifold)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      integer manifold
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj
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
                     qf = qfine(i2(m),j2(m),mq)
                     kf = areafine(i2(m),j2(m))
                     sum = sum + kf*qf
                  enddo
                  kc = areacoarse(i1,j1)
                  qcoarse(i1,j1,mq) = sum/kc
               else
                  sum = 0
                  do m = 0,r2-1
                     qf = qfine(i2(m),j2(m),mq)
                     sum = sum + qf
                  enddo
                  qcoarse(i1,j1,mq) = sum/r2
               endif
            enddo
         enddo
      enddo
      end


