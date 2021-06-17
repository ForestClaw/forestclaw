c     # ----------------------------------------------------------------------
c     # This handles the boundary conditions at the block
c     # corners for the pillow sphere.
c     # ----------------------------------------------------------------------

      subroutine fclaw2d_pillow46_copy_block_corner
     &      (mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(1-ibc,1-jbc,mq) = qneighbor(ibc,jbc,mq)
                  qneighbor(1-ibc,1-jbc,mq) = qthis(ibc,jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mx+ibc,1-jbc,mq) = qneighbor(mx+1-ibc,jbc,mq)
                  qneighbor(mx+ibc,1-jbc,mq) = qthis(mx+1-ibc,jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(1-ibc,my+jbc,mq) = qneighbor(ibc,my+1-jbc,mq)
                  qneighbor(1-ibc,my+jbc,mq) = qthis(ibc,my+1-jbc,mq)
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(mx+1-ibc,my+1-jbc,mq)
                  qneighbor(mx+ibc,my+jbc,mq) =
     &                  qthis(mx+1-ibc,my+1-jbc,mq)
               enddo
            enddo
         endif
      enddo
      end

      subroutine fclaw2d_pillow46_average_block_corner
     &      (mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, areacoarse, areafine,
     &      icorner,iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock, refratio
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc, ii, jj, ifine, jfine
      double precision sum,qf,kf,kc


      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(ifine,jfine,mq)
                        kf = areafine(ifine,jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(1-ibc,1-ibc)
                  qcoarse(1-ibc,1-jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(mx+1-ifine,jfine,mq)
                        kf = areafine(mx+1-ifine,jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(mx+ibc,1-jbc)
                  qcoarse(mx+ibc,1-jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(ifine,my+1-jfine,mq)
                        kf = areafine(ifine,my+1-jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(1-ibc,my+jbc)
                  qcoarse(1-ibc,my+jbc,mq) = sum/kc
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        qf = qfine(mx+1-ifine,my+1-jfine,mq)
                        kf = areafine(mx+1-ifine,my+1-jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(mx+ibc,my+jbc)
                  qcoarse(mx+ibc,my+jbc,mq) = sum/kc
               enddo
            enddo
         endif
      enddo

      end


      subroutine fclaw2d_pillow46_interpolate_block_corner
     &      (mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, icorner_coarse, iblock)
      implicit none
      integer mx, my, mbc, meqn, icorner_coarse, iblock, refratio
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc
      integer ic, jc, mth
      double precision gradx, grady, shiftx, shifty
      double precision sl, sr, qc, value
      double precision fclaw2d_clawpatch_compute_slopes

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

c     # This may not even matter
      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)
         sl = (qc - qcoarse(ic-1,jc,mq))
         sr = (qcoarse(ic+1,jc,mq) - qc)
         gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(ic,jc-1,mq))
         sr = (qcoarse(ic,jc+1,mq) - qc)
         grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

c        # Loop over fine grid ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = qc + shiftx*gradx + shifty*grady
               if (icorner_coarse .eq. 0) then
                  qfine(1-ibc,1-jbc,mq) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(mx+ibc,1-jbc,mq) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(1-ibc,my+jbc,mq) = value
               elseif (icorner_coarse .eq. 3) then
                  qfine(mx+ibc,my+jbc,mq) = value
               endif
            enddo
         enddo
      enddo
      end
