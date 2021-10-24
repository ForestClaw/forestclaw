c ----------------------------------------------------------------------
c> @file
c> This handles the boundary conditions at the block
c> corners for the pillow sphere.
c ----------------------------------------------------------------------

c--------------------------------------------------------------------
c> @brief @copybrief ::pillow_fort_copy_block_corner_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::pillow_fort_copy_block_corner_t
c--------------------------------------------------------------------
      subroutine fclaw2d_pillow5_copy_block_corner
     &      (mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock
      double precision qthis(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qneighbor(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq, ibc, jbc

      do mq = 1,meqn
         if (icorner .eq. 0) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mq,1-ibc,1-jbc) = qneighbor(mq,ibc,jbc)
                  qneighbor(mq,1-ibc,1-jbc) = qthis(mq,ibc,jbc)
               enddo
            enddo
         elseif (icorner .eq. 1) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mq,mx+ibc,1-jbc) = qneighbor(mq,mx+1-ibc,jbc)
                  qneighbor(mq,mx+ibc,1-jbc) = qthis(mq,mx+1-ibc,jbc)
               enddo
            enddo
         elseif (icorner .eq. 2) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mq,1-ibc,my+jbc) = qneighbor(mq,ibc,my+1-jbc)
                  qneighbor(mq,1-ibc,my+jbc) = qthis(mq,ibc,my+1-jbc)
               enddo
            enddo
         elseif (icorner .eq. 3) then
            do ibc = 1,mbc
               do jbc = 1,mbc
                  qthis(mq,mx+ibc,my+jbc) =
     &                  qneighbor(mq,mx+1-ibc,my+1-jbc)
                  qneighbor(mq,mx+ibc,my+jbc) =
     &                  qthis(mq,mx+1-ibc,my+1-jbc)
               enddo
            enddo
         endif
      enddo
      end

c--------------------------------------------------------------------
c> @brief @copybrief ::pillow_fort_average_block_corner_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::pillow_fort_average_block_corner_t
c--------------------------------------------------------------------
      subroutine fclaw2d_pillow5_average_block_corner
     &      (mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, areacoarse, areafine,
     &      icorner,iblock)
      implicit none

      integer mx, my, mbc, meqn, icorner, iblock, refratio
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
                        qf = qfine(mq,ifine,jfine)
                        kf = areafine(ifine,jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(1-ibc,1-ibc)
                  qcoarse(mq,1-ibc,1-jbc) = sum/kc
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
                        qf = qfine(mq,mx+1-ifine,jfine)
                        kf = areafine(mx+1-ifine,jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(mx+ibc,1-jbc)
                  qcoarse(mq,mx+ibc,1-jbc) = sum/kc
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
                        qf = qfine(mq,ifine,my+1-jfine)
                        kf = areafine(ifine,my+1-jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(1-ibc,my+jbc)
                  qcoarse(mq,1-ibc,my+jbc) = sum/kc
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
                        qf = qfine(mq,mx+1-ifine,my+1-jfine)
                        kf = areafine(mx+1-ifine,my+1-jfine)
                        sum = sum + kf*qf
                     enddo
                  enddo
                  kc = areacoarse(mx+ibc,my+jbc)
                  qcoarse(mq,mx+ibc,my+jbc) = sum/kc
               enddo
            enddo
         endif
      enddo

      end


c--------------------------------------------------------------------
c> @brief @copybrief ::pillow_fort_average_block_corner_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::pillow_fort_average_block_corner_t
c--------------------------------------------------------------------
      subroutine fclaw2d_pillow5_interpolate_block_corner
     &      (mx,my,mbc,meqn,
     &      refratio, qcoarse, qfine, icorner_coarse, iblock)
      implicit none
      integer mx, my, mbc, meqn, icorner_coarse, iblock, refratio
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

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
      else
         write(6,*) 'pillow : icorner_coarse has unexpected value'
         write(6,*) 'icorner_coarse : ', icorner_coarse
         stop
      endif

c     # This may not even matter
      do mq = 1,meqn
         qc = qcoarse(mq,ic,jc)
         sl = (qc - qcoarse(mq,ic-1,jc))
         sr = (qcoarse(mq,ic+1,jc) - qc)
         gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

         sl = (qc - qcoarse(mq,ic,jc-1))
         sr = (qcoarse(mq,ic,jc+1) - qc)
         grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

c        # Loop over fine grid ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = qc + shiftx*gradx + shifty*grady
               if (icorner_coarse .eq. 0) then
                  qfine(mq,1-ibc,1-jbc) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(mq,mx+ibc,1-jbc) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(mq,1-ibc,my+jbc) = value
               elseif (icorner_coarse .eq. 3) then
                  qfine(mq,mx+ibc,my+jbc) = value
               endif
            enddo
         enddo
      enddo
      end
