      subroutine fc2d_geoclaw_fort_interpolate2fine(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mbathy,
     &      p4est_refineFactor,refratio,igrid)

      use geoclaw_module, ONLY:dry_tolerance, sea_level
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
      double precision etabarc(-1:1, -1:1), h, b, u, hfsum
      double precision coarseumin, coarseumax
      logical redefine

      hfsum = 0.d0
c     # Use minmod to calculate slope.
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
               ! Calculate surface elevation eta using dry limiting
               if (mq .eq. 1) then
                  do ii = -1, 1
                     do jj = -1, 1
                        h = qcoarse(1,ic+ii,jc+jj)
                        b = aux_coarse(mbathy,ic+ii,jc+jj)
                        if (h < dry_tolerance) then
                           etabarc(ii,jj) = sea_level
                        else
                           etabarc(ii,jj) = h + b
                        endif
                     enddo
                  enddo
c                 # Interpolate sea surface height rather than just the
c                 # water column height.
                  qc = etabarc(0,0)

                  sl = qc - etabarc(-1,0)
                  sr = etabarc(1,0) - qc
                  gradx = compute_slopes(sl,sr,mth)

                  sl = qc - etabarc(0,-1)
                  sr = etabarc(0,1) - qc
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
                        qfine(mq,iff,jf) = max(qf,0.0)
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
                        hfsum = hfsum + qfine(mq,iff,jf)
                     enddo
                  enddo


c------------check to make sure we are not creating any new extrema

c                 # calculate coarse cells' velocity
c                  do ii = -1,1
c                     do jj = -1,1
c                        if (qcoarse(1,ic+ii,jc+jj) .eq. 0.d0) then
c                           uc(ii,jj) = 0.d0
c                        else
c                           uc(ii,jj) = qcoarse(mq,ic+ii,jc+jj)/
c     &                                 qcoarse(1,ic+ii,jc+jj)
c                        endif
c                     enddo
c                  enddo
c                 # find the maximum/minimum velocities among coarse cells
c                  coarseumax = -1d99
c                  coarseumin = 1d99
c                  do ii = -1,1
c                     do jj = -1,1
c                        coarseumax = max(coarseumax,uc(ii,jj))
c                        coarseumin = min(coarseumin,uc(ii,jj))
c                     enddo
c                  enddo

c                  redefine = .false.
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        iff = (i-1)*refratio + ii
c                        jf = (j-1)*refratio + jj
c                        if (qfine(1,iff,jf) .eq. 0.d0) then
c                           uf = 0.d0
c                        else
c                           uf = qfine(mq,iff,jf)/qfine(1,iff,jf)
c                        endif
c                        if (uf .gt. coarseumax .or. uf .lt. coarseumin)
c     &                        then
c                           redefine = .true.
c                        endif
c                     enddo
c                  enddo

c                  if (redefine) then
c                     u = qcoarse(mq,ic,jc)/
c     &                 max(qcoarse(1,ic,jc), hfsum/(refratio*refratio))
c                     do ii = 1,refratio
c                        do jj = 1,refratio
c                           iff = (i-1)*refratio + ii
c                           jf = (j-1)*refratio + jj
c                           qfine(mq,iff,jf) = qfine(1,iff,jf)*u
c                        enddo
c                     enddo
c                  endif
c------end of checking to make sure we are not creating any new extrema

               endif
c------end of interpolation
            enddo
         enddo
      enddo

      end
