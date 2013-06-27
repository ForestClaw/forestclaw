      subroutine interpolate_geo(mx,my,mbc,meqn, qcoarse,
     &      qfine, maux,aux_coarse,aux_fine,
     &      p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid, maux
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision aux_coarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision aux_fine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

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
                  qc = qcoarse(ic,jc,mq) + aux_coarse(ic,jc,19)

                  sl = (qc - qcoarse(ic-1,jc,mq) -
     &                  aux_coarse(ic-1,jc,19))
                  sr = (qcoarse(ic+1,jc,mq)+aux_coarse(ic+1,jc,19) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(ic,jc-1,mq) -
     &                  aux_coarse(ic,jc-1,19))
                  sr = (qcoarse(ic,jc+1,mq)+aux_coarse(ic,jc+1,19) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # Fill in fine grid values from coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qf = qc + shiftx*gradx + shifty*grady -
     &                        aux_fine(iff,jf,19)
                        qfine(iff,jf,mq) = qf
                     enddo
                  enddo
               else
c                 # interpolate momentum components in the usual way.
c                 # But then make sure that no new extrema are created.
                  qc = qcoarse(ic,jc,mq)

                  sl = (qc - qcoarse(ic-1,jc,mq))
                  sr = (qcoarse(ic+1,jc,mq) - qc)
                  gradx = compute_slopes(sl,sr,mth)

                  sl = (qc - qcoarse(ic,jc-1,mq))
                  sr = (qcoarse(ic,jc+1,mq) - qc)
                  grady = compute_slopes(sl,sr,mth)

c                 # Fill in refined values on coarse grid cell (ic,jc)
                  do ii = 1,refratio
                     do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        qfine(iff,jf,mq) =
     &                        qc + shiftx*gradx + shifty*grady
                     enddo
                  enddo

c                 # check to make sure we are not creating any new extrema
                  do ii = -1,1
                     do jj = -1,1
                        uc(ii,jj) = qcoarse(ic+ii,jc+jj,mq)/
     &                        qcoarse(ic+ii,jc+jj,1)
                     enddo
                  enddo

                  coarseumax = -1d99
                  coarseumin = 1d99
                  do ii = -1,1
                     do jj = -1,1
                        coarseumax = max(coarseumax,uc(ii,jj))
                        coarseumin = min(coarseumin,uc(ii,jj))
                     enddo
                  enddo

                  redefine = .false.
                  do ii = 1,refratio
                     do jj = 1,refratio
                        iff = (i-1)*refratio + ii
                        jf = (j-1)*refratio + jj
                        uf = qfine(iff,jf,mq)/qfine(iff,jf,1)
                        if (uf .gt. coarseumax .or. uf .lt. coarseumin)
     &                        then
                           redefine = .true.
                        endif
                     enddo
                  enddo

                  if (redefine) then
                     do ii = 1,refratio
                        do jj = 1,refratio
                           iff = (i-1)*refratio + ii
                           jf = (j-1)*refratio + jj
                           qfine(iff,jf,mq) = qfine(iff,jf,1)*uc(0,0)
                        enddo
                     enddo
                  endif

               endif
            enddo
         enddo
      enddo

      end
