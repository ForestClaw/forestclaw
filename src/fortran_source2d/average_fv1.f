c     # ----------------------------------------------------------
c     # Averaging routines - (i,j,mq) ordering
c     # ----------------------------------------------------------
c     # average_face_ghost
c     # average_corner_ghost
c     # average_to_coarse_patch
c     # ------------------------------------------------------------------

c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iface_coarse'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost(mx,my,mbc,meqn,
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

      integer mq,r2, m, m1
      integer ic, ic1, ibc
      integer jc, jc1, jbc

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      double precision kc

      logical is_valid_average, skip_this_grid

      is_manifold = manifold .eq. 1

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

c     # Average fine grid onto coarse grid
      do mq = 1,meqn
         if (idir .eq. 0) then
            do jc = 1,my
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer

                  if (iface_coarse .eq. 0) then
                     ic = 1-ibc
                  elseif (iface_coarse .eq. 1) then
                     ic = mx+ibc
                  endif

c                 # New code
                  ic1 = ic
                  jc1 = jc
                  call fclaw2d_transform_face_half(ic1,jc1,i2,j2,
     &                  transform_cptr)
                  skip_this_grid = .false.
                  do m = 0,r2-1
                     if (.not. is_valid_average(i2(m),j2(m),mx,my))
     &                     then
                        skip_this_grid = .true.
                     endif
                  enddo

c                  if (iface_coarse .eq. 0 .and. ic .eq. 0 .and.
c     &                  igrid .eq. 0) then
c                     write(6,'(A,2I4)') 'ic,jc : ', ic,jc
c                     write(6,'(A,4I4)') 'i2    : ', (i2(m),m=0,3)
c                     write(6,'(A,4I4)') 'j2    : ', (j2(m),m=0,3)
c                     write(6,*) ' '
c                  endif
                  if (.not. skip_this_grid) then
                     if (is_manifold) then
                        sum = 0
                        do m = 0,r2-1
                           qf = qfine(i2(m),j2(m),mq)
                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                        enddo
                        kc = areacoarse(ic,jc)
                        qcoarse(ic,jc,mq) = sum/kc
                     else
                        sum = 0
                        do m = 0,r2-1
                           sum = sum + qfine(i2(m),j2(m),mq)
                        enddo
                        qcoarse(ic,jc,mq) = sum/r2
                     endif
                  endif
               enddo
            enddo
         else
c           # idir = 1 (faces 2,3)
            do ic = 1,mx
               do jbc = 1,mbc

                  if (iface_coarse .eq. 2) then
                     jc = 1-jbc
                  elseif (iface_coarse .eq. 3) then
                     jc = my+jbc
                  endif

                  ic1 = ic
                  jc1 = jc
                  call fclaw2d_transform_face_half(ic1,jc1,i2,j2,
     &                  transform_cptr)
                  skip_this_grid = .false.
                  do m = 0,r2-1
                     if (.not. is_valid_average(i2(m),j2(m),mx,my))
     &                     then
                        skip_this_grid = .true.
                     endif
                  enddo
                  if (.not. skip_this_grid) then
                     if (is_manifold) then
                        sum = 0
                        do m = 0,r2-1
                           qf = qfine(i2(m),j2(m),mq)
                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                        enddo
                        kc = areacoarse(ic,jc)
                        qcoarse(ic,jc,mq) = sum/kc
                     else
                        sum = 0
                        do m = 0,r2-1
                           sum = sum + qfine(i2(m),j2(m),mq)
                        enddo
                        qcoarse(ic,jc,mq) = sum/r2
                     endif
                  endif
               enddo
            enddo
         endif
      enddo

      end

      logical function is_valid_average(i,j,mx,my)
      implicit none

      integer i,j,mx,my
      logical i1, j1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my

      is_valid_average = i1 .and. j1

      end


c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_corner_ghost(mx,my,mbc,meqn,
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

      integer i,j,ibc,jbc,ii,jj,mq,r2
      integer ifine, jfine
      logical is_manifold
      double precision qf,kf, kc

c     # This should be refratio*refratio.
      integer i1,j1,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

      is_manifold = manifold .eq. 1

      r2 = refratio*refratio
      do mq = 1,meqn
c        # Loop over four corner cells on coarse grid
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Average fine grid corners onto coarse grid ghost corners
               if (icorner_coarse .eq. 0) then
    1             i1 = 1-ibc
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

c              # Again, a fake routine until the real one is
c              # available (be sure to pass in (i1,j1)
               call fclaw2d_transform_corner_half(i1,j1,i2,j2,
     &               transform_cptr)
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
                     sum = sum + qfine(i2(m),j2(m),mq)
                  enddo
                  qcoarse(i1,j1,mq) = sum/r2
               endif


            enddo
         enddo
      enddo


      end

      subroutine average_to_coarse_patch(mx,my,mbc,meqn,
     &      qcoarse,qfine,
     &      areacoarse, areafine,
     &      p4est_refineFactor,
     &      refratio, igrid,manifold)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      integer manifold
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

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
