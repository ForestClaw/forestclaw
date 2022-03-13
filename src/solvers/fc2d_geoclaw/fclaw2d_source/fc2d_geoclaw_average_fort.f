c     > \file
c     > \defgroup Averaging Average fine grids to a coarse grid
c     > Average cells from coarse grid to fine grid.
c     >
c     > Routines described here average are used to fill coarse grid ghost
c     > cells, and average sibling grids onto a parent grid.  Indices
c     > for cells at block boundaries are transformed using encodings
c     > stored in `transform_cptr`.
c     >
c     > \param [in] mx,my       Number of cells in x,y direction
c     > \param [in] mbc      Number of ghost cells
c     > \param [in] meqn     Number of equations
c     > \param [in] qcoarse,qfine  Solution on coarse,fine grid
c     > \param [in] areacoarse,areafine  Area of mesh cells on coarse,fine
c     grids.
c     > \param [in] idir     Face orientation - 0 for x-faces; 1 for y-faces
c     [0-1]
c     > \param [in] iface    Face number of fine grid [0-3].
c     > \param [in] iface_coarse Face number of coarse grid [0-3].
c     > \param [in] num_neighbors Number of fine grid neighbors [2].
c     > \param [in] refratio  Refinement ratio between coarse and fine grids
c     [2].
c     > \param [in] manifold  Flag indicating whether we are on mapped grid [0
c     -1].
c     > \param [in] transform_cptr  Encoding for indices at block boundaries
c     (C only).

c     > \ingroup Averaging
c     > Average fine ghost cell values.
c     >
c     > Average fine grid interior values to neighboring ghost cell values of
c     > the coarse grid.
      subroutine fc2d_geoclaw_fort_average_face(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mcapa,
     &      idir,iface_coarse,igrid, transform_cptr)

      use geoclaw_module, only: dry_tolerance

      implicit none

      integer mx,my,mbc,meqn,igrid,idir,iface_coarse
      integer mcapa,maux
      integer*8 transform_cptr
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

c     # these will be empty if we are not on a manifold.
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision sum, qf, kf
      logical is_manifold

      integer mq,r2, m !, m1
      integer ic,  ibc !ic1,
      integer jc,  jbc !jc1,

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      !double precision kc

      logical fclaw2d_clawpatch_is_valid_average, skip_this_grid
      double precision af_sum, qv(0:rr2-1)

      double precision etasum, hsum, husum, hvsum, etaav, hav
      double precision hc, huc, hvc
      double precision hf, huf, hvf, bf, etaf
      double precision capac, capa
      integer nwet

      integer mbathy, refratio, num_neighbors

      mbathy = 1
      refratio = 2
      num_neighbors = 2

c     # Even though Geoclaw is on a manifold, we don't handle it 
c     # in the usual forestclaw way      
      is_manifold = .false.

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

c     # Average fine grid onto coarse grid
      if (idir .eq. 0) then
         do jc = 1,my
            do ibc = 1,mbc
c              # ibc = 1 corresponds to first layer of ghost cells, and
c              # ibc = 2 corresponds to the second layer

               if (iface_coarse .eq. 0) then
                  ic = 1-ibc
               elseif (iface_coarse .eq. 1) then
                  ic = mx+ibc
               endif

               call fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &               transform_cptr)
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
                  if (.not. 
     &        fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo

               if (.not. skip_this_grid) then
                  if (is_manifold) then
                     write(6,'(A,A,A)') 'fc2d_geoclaw_average_fort : ',
     &                   'We should not be here;  manifold case is ',
     &                   'not handled here in geoclaw'
                     stop
                     do mq = 1,meqn
                        sum = 0
                        af_sum = 0
                        do m = 0,r2-1
                           qf = qfine(mq,i2(m),j2(m))
                           qv(m) = qf
c                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                           af_sum = af_sum + kf
                        enddo
c                       # ----------------------------------------
c                       # At block seams, the coarse grid mesh cell
c                       # areas may not have been computed using
c                       # the correct metrics.
c                       # ----------------------------------------
c                       kc = areacoarse(ic,jc)
c                       qcoarse(mq,ic,jc) = sum/kc

c                       # Use areas of the fine grid mesh cells instead.
                        qcoarse(mq,ic,jc) = sum/af_sum
                     enddo
                  else
                     if (mcapa .eq. 0) then
                        capac=1.0d0
                     else
                        capac=aux_coarse(mcapa,ic,jc)
                     endif 
                     etasum = 0.d0
                     hsum   = 0.d0
                     husum  = 0.d0
                     hvsum  = 0.d0

                     nwet   = 0 
c                    sum = 0
                     do m = 0,r2-1
c                       sum = sum + qfine(mq,i2(m),j2(m))
c                       qcoarse(mq,ic,jc) = sum/dble(r2)
                        if (mcapa .eq. 0) then
                           capa=1.0d0
                        else
                           capa=aux_fine(mcapa,i2(m),j2(m))
                        endif
                        hf = qfine(1,i2(m),j2(m))*capa
                        bf = aux_fine(mbathy,i2(m),j2(m))*capa
                        huf= qfine(2,i2(m),j2(m))*capa
                        hvf= qfine(3,i2(m),j2(m))*capa
                        if (hf > dry_tolerance) then
                           etaf = hf+bf
                           nwet=nwet+1
                        else
                           etaf = 0.d0
                           huf=0.d0
                           hvf=0.d0
                        endif
                        hsum   = hsum + hf
                        husum  = husum + huf
                        hvsum  = hvsum + hvf
                        etasum = etasum + etaf
                     enddo
                     if (nwet.gt.0) then
                        etaav=etasum/dble(nwet)
                        hav=hsum/dble(nwet)
c                       hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
                        hc=min(hav,(max(etaav-
     &                        aux_coarse(mbathy,ic,jc)*capac,0.d0)))
c                       huc=(min(hav,hc)/hsum)*husum
c                       hvc=(min(hav,hc)/hsum)*hvsum
                        huc=(hc/hsum)*husum
                        hvc=(hc/hsum)*hvsum
                     else
                        hc=0.d0
                        huc=0.d0
                        hvc=0.d0
                     endif
                     qcoarse(1,ic,jc) = hc / capac
                     qcoarse(2,ic,jc) = huc / capac
                     qcoarse(3,ic,jc) = hvc / capac
                  endif
               endif
            enddo
         enddo
      else
c        # idir = 1 (faces 2,3)
         do jbc = 1,mbc
            do ic = 1,mx

               if (iface_coarse .eq. 2) then
                  jc = 1-jbc
               elseif (iface_coarse .eq. 3) then
                  jc = my+jbc
               endif

               call fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &               transform_cptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. 
     &         fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my))
     &                  then
                     skip_this_grid = .true.
                  endif
               enddo
               if (.not. skip_this_grid) then
                  if (is_manifold) then
                     do mq = 1,meqn
                        sum = 0
                        af_sum = 0
                        do m = 0,r2-1
                           qf = qfine(mq,i2(m),j2(m))
c                           kf = areafine(i2(m),j2(m))
                           sum = sum + qf*kf
                           af_sum = af_sum + kf
                        enddo
c                        kc = areacoarse(ic,jc)
c                       qcoarse(mq,ic,jc) = sum/kc
                        qcoarse(mq,ic,jc) = sum/af_sum
                     enddo
                  else
                     if (mcapa .eq. 0) then
                        capac=1.0d0
                     else
                        capac=aux_coarse(mcapa,ic,jc)
                     endif 
                     etasum = 0.d0
                     hsum   = 0.d0
                     husum  = 0.d0
                     hvsum  = 0.d0

                     nwet   = 0                   
                     do m = 0,r2-1
c                       sum = sum + qfine(mq,i2(m),j2(m))
c                       qcoarse(mq,ic,jc) = sum/dble(r2)
                        if (mcapa .eq. 0) then
                           capa=1.0d0
                        else
                           capa=aux_fine(mcapa,i2(m),j2(m))
                        endif
                        hf = qfine(1,i2(m),j2(m))*capa
                        bf = aux_fine(mbathy,i2(m),j2(m))*capa
                        huf= qfine(2,i2(m),j2(m))*capa
                        hvf= qfine(3,i2(m),j2(m))*capa
                        if (hf > dry_tolerance) then
                           etaf = hf+bf
                           nwet=nwet+1
                        else
                           etaf = 0.d0
                           huf=0.d0
                           hvf=0.d0
                        endif
                        hsum   = hsum + hf
                        husum  = husum + huf
                        hvsum  = hvsum + hvf
                        etasum = etasum + etaf
                     enddo
                     if (nwet.gt.0) then
                        etaav=etasum/dble(nwet)
                        hav=hsum/dble(nwet)
c                       hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
     &
                        hc=min(hav,(max(etaav-
     &                        aux_coarse(mbathy,ic,jc)*capac,0.d0)))
c                       huc=(min(hav,hc)/hsum)*husum
c                       hvc=(min(hav,hc)/hsum)*hvsum
                        huc=(hc/hsum)*husum
                        hvc=(hc/hsum)*hvsum
                     else
                        hc=0.d0
                        huc=0.d0
                        hvc=0.d0
                     endif
                     qcoarse(1,ic,jc) = hc / capac
                     qcoarse(2,ic,jc) = huc / capac
                     qcoarse(3,ic,jc) = hvc / capac


c                     sum = 0
c                     do m = 0,r2-1
c                        sum = sum + qfine(mq,i2(m),j2(m))
c                     enddo
c                     qcoarse(mq,ic,jc) = sum/dble(r2)
                  endif                 !! manifold loop
               endif                    !! skip grid loop
            enddo
         enddo
      endif

      end



c> \ingroup Averaging
c> Average across corners.
      subroutine fc2d_geoclaw_fort_average_corner(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mcapa,
     &      icorner_coarse,transform_cptr)
      
      use geoclaw_module, only: dry_tolerance
      
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse
      integer*8 transform_cptr
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

c     # these will be empty if we are not on a manifold.
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision sum

      integer ibc,jbc,mq,r2 !i,j,jj,ii,
      !integer  jfine !ifine
      logical is_manifold
      double precision qf,kf !, kc

c     # This should be refratio*refratio.
      integer i1,j1,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      double precision af_sum

      double precision etasum, hsum, husum, hvsum, etaav, hav
      double precision hc, huc, hvc
      double precision hf, huf, hvf, bf, etaf
      double precision capac, capa
      integer nwet

      integer maux,mcapa,mbathy, refratio

      refratio = 2
      mbathy = 1

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

c     # Even though Geoclaw is on a manifold, we don't handle it 
c     # in the usual forestclaw way      
      is_manifold = .false.

      r2 = refratio*refratio
c     # Loop over four corner cells on coarse grid
      do ibc = 1,mbc
         do jbc = 1,mbc
c           # Average fine grid corners onto coarse grid ghost corners
            if (icorner_coarse .eq. 0) then
    1          i1 = 1-ibc
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
                     qf = qfine(mq,i2(m),j2(m))
c                     kf = areafine(i2(m),j2(m))
                     sum = sum + kf*qf
                     af_sum = af_sum + kf
                  enddo
c                  kc = areacoarse(i1,j1)
c                 qcoarse(mq,i1,j1) = sum/kc
                  qcoarse(mq,i1,j1) = sum/af_sum
               enddo
            else
               if (mcapa .eq. 0) then
                  capac = 1.0d0
               else
                  capac = aux_coarse(mcapa,i1,j1)
               endif 
               etasum = 0.d0
               hsum   = 0.d0
               husum  = 0.d0
               hvsum  = 0.d0

               nwet   = 0                   
               do m = 0,r2-1
c                       sum = sum + qfine(mq,i2(m),j2(m))
c                       qcoarse(mq,ic,jc) = sum/dble(r2)
                  if (mcapa .eq. 0) then
                     capa = 1.0d0
                  else
                     capa = aux_fine(mcapa,i2(m),j2(m))
                  endif
                  hf = qfine(1,i2(m),j2(m))*capa
                  bf = aux_fine(mbathy,i2(m),j2(m))*capa
                  huf= qfine(2,i2(m),j2(m))*capa
                  hvf= qfine(3,i2(m),j2(m))*capa
                  if (hf > dry_tolerance) then
                     etaf = hf+bf
                     nwet = nwet+1
                  else
                     etaf = 0.d0
                     huf = 0.d0
                     hvf = 0.d0
                  endif
                  hsum   = hsum + hf
                  husum  = husum + huf
                  hvsum  = hvsum + hvf
                  etasum = etasum + etaf
               enddo
               if (nwet.gt.0) then
                  etaav = etasum/dble(nwet)
                  hav = hsum/dble(nwet)
c                       hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
                  hc = min(hav,(max(etaav-
     &                        aux_coarse(mbathy,i1,j1)*capac,0.d0)))
c                       huc=(min(hav,hc)/hsum)*husum
c                       hvc=(min(hav,hc)/hsum)*hvsum
                  huc = (hc/hsum)*husum
                  hvc = (hc/hsum)*hvsum
               else
                  hc = 0.d0
                  huc = 0.d0
                  hvc = 0.d0
               endif
               qcoarse(1,i1,j1) = hc / capac
               qcoarse(2,i1,j1) = huc / capac
               qcoarse(3,i1,j1) = hvc / capac
c               do mq = 1,meqn
c                  sum = 0
c                  do m = 0,r2-1
c                     sum = sum + qfine(mq,i2(m),j2(m))
c                  enddo
c                  qcoarse(mq,i1,j1) = sum/dble(r2)
c               enddo
            endif
         enddo
      enddo

      end

      subroutine fc2d_geoclaw_fort_average2coarse(mx,my,mbc,meqn,
     &      qcoarse,qfine,maux,aux_coarse,aux_fine,mcapa,igrid)

      use geoclaw_module, only: dry_tolerance      
      
      implicit none

      integer mx,my,mbc,meqn
      integer igrid,maux,mcapa,nwet 
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj !, ifine, jfine
      double precision etasum, hsum, husum, hvsum, etaav, hav
      double precision hc, huc, hvc
      double precision hf, huf, hvf, bf, etaf
      double precision capac, capa

      integer refratio, p4est_refineFactor,mbathy

c     # This should be refratio*refratio.
      integer i1,j1,r2,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      refratio = 2
      p4est_refineFactor = 2
      mbathy = 1

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
      jc_add = jg*my/p4est_refineFactor

      r2 = refratio * refratio

c     # First loop over quadrant (i1,i2)x(j1,j2) of the coarse grid
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

            if (mcapa .eq. 0) then
               capac=1.0d0
            else
               capac=aux_coarse(mcapa,i1,j1)
            endif  

            etasum = 0.d0
            hsum   = 0.d0
            husum  = 0.d0
            hvsum  = 0.d0

            nwet   = 0
c           # loop over the fine grids
            do m = 0,r2-1
               if (mcapa .eq. 0) then
                  capa=1.0d0
               else
                  capa=aux_fine(mcapa,i2(m),j2(m))
                  endif
               hf = qfine(1,i2(m),j2(m))*capa
               bf = aux_fine(mbathy,i2(m),j2(m))*capa
               huf= qfine(2,i2(m),j2(m))*capa 
               hvf= qfine(3,i2(m),j2(m))*capa
               if (hf > dry_tolerance) then
                  etaf = hf+bf
                  nwet=nwet+1
               else
                  etaf = 0.d0
                  huf=0.d0
                  hvf=0.d0
                  endif
               hsum   = hsum + hf
               husum  = husum + huf
               hvsum  = hvsum + hvf
               etasum = etasum + etaf 
               enddo          
            if (nwet.gt.0) then
               etaav=etasum/dble(nwet)
               hav=hsum/dble(nwet)
c              hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
               hc=min(hav,(max(etaav-
     &             aux_coarse(1,i1,j1)*capac,0.d0)))
c               huc=(min(hav,hc)/hsum)*husum
c               hvc=(min(hav,hc)/hsum)*hvsum
               huc=(hc/hsum)*husum
               hvc=(hc/hsum)*hvsum
            else
               hc=0.d0
               huc=0.d0
               hvc=0.d0
            endif      
            qcoarse(1,i1,j1) = hc / capac 
            qcoarse(2,i1,j1) = huc / capac 
            qcoarse(3,i1,j1) = hvc / capac 
            enddo
         enddo

      end

