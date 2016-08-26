      subroutine geoclaw_average2coarse(mx,my,mbc,meqn,qcoarse,
     &      qfine,maux,aux_coarse,aux_fine,mcapa,mbathy,
     &      p4est_refineFactor,refratio,igrid)

      use geoclaw_module, only: dry_tolerance      
      
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid,maux,mbathy,mcapa,nwet 
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      double precision etasum, hsum, husum, hvsum, etaav, hav
      double precision hc, huc, hvc
      double precision hf, huf, hvf, bf, etaf
      double precision capac, capa

c     # This should be refratio*refratio.
      integer i1,j1,r2,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

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
     &             aux_coarse(mbathy,i1,j1)*capac,0.d0)))
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
