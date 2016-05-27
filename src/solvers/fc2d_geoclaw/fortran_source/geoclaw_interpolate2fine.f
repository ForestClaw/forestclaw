      subroutine geoclaw_interpolate2fine(mx,my,mbc,meqn,qcoarse,
     &      qfine,maux,aux_coarse,aux_fine,mbathy,
     &      p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid, maux,mbathy
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      integer ic,jc,ic_add, jc_add, iff, jf
      double precision qc, qf, shiftx, shifty, sl, sr, gradx, grady
      double precision compute_slopes, uc(-1:1,-1:1), uf
      double precision coarseumin, coarseumax
      logical redefine


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
c        # First loop over quadrant (i1,i2)x(j1,j2) of the coarse grid
         do i = i1,i2
            do j = j1,j2
               ic = i + ic_add
               jc = j + jc_add
               if (mq .eq. 1) then
c                 # Interpolate sea surface height rather than just the
c                 # water column height.
                  qc = qcoarse(mq,ic,jc) + aux_coarse(mbathy,ic,jc)

                  sl = (qc - qcoarse(mq,ic-1,jc) -
     &                  aux_coarse(mbathy,ic-1,jc))
                  sr = (qcoarse(mq,ic+1,jc) +
     &                  aux_coarse(mbathy,ic+1,jc) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(mq,ic,jc-1) -
     &                  aux_coarse(mbathy,ic,jc-1))
                  sr = (qcoarse(mq,ic,jc+1)+
     &                  aux_coarse(mbathy,ic,jc+1) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # Fill in fine grid values from coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qf = qc + shiftx*gradx + shifty*grady -
     &                        aux_fine(mbathy,iff,jf)
                        qfine(mq,iff,jf) = qf
                     enddo
                  enddo
               else
c                 # interpolate momentum components in the usual way.
c                 # But then make sure that no new extrema are created.
                  qc = qcoarse(mq,ic,jc)

                  sl = (qc - qcoarse(mq,ic-1,jc))
                  sr = (qcoarse(mq,ic+1,jc) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(mq,ic,jc-1))
                  sr = (qcoarse(mq,ic,jc+1) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # Fill in refined values on coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qfine(mq,iff,jf) =
     &                        qc + shiftx*gradx + shifty*grady
                     enddo
                  enddo

cc                 # check to make sure we are not creating any new extrema
c                  do ii = -1,1
c                     do jj = -1,1
c                        uc(ii,jj) = qcoarse(mq,ic+ii,jc+jj)/
c     &                        qcoarse(1,ic+ii,jc+jj)
c                     enddo
c                  enddo
c
c                  coarseumax = -1d99
c                  coarseumin = 1d99
c                  do ii = -1,1
c                     do jj = -1,1
c                        coarseumax = max(coarseumax,uc(ii,jj))
c                        coarseumin = min(coarseumin,uc(ii,jj))
c                     enddo
c                  enddo
c
c                  redefine = .false.
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        iff = (i-1)*refratio + ii
c                        jf = (j-1)*refratio + jj
c                        uf = qfine(mq,iff,jf)/qfine(1,iff,jf)
c                        if (uf .gt. coarseumax .or. uf .lt. coarseumin)
c     &                        then
c                           redefine = .true.
c                        endif
c                     enddo
c                  enddo
c
c                  if (redefine) then
c                     do ii = 1,refratio
c                        do jj = 1,refratio
c                           iff = (i-1)*refratio + ii
c                           jf = (j-1)*refratio + jj
c                           qfine(mq,iff,jf) = qfine(1,iff,jf)*uc(0,0)
c                        enddo
c                     enddo
c                  endif

               endif
            enddo
         enddo
      enddo

      end
