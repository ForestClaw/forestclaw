c    # ----------------------------------------------------------------
c    # This file contains routines which accumulate fluxes, waves 
c    # and add corrections to coarse grid at edges of both same size
c    # grids and coarse/fine interfaces. 
c    # 
c    # 1. clawpack46_time_sync_store_flux (before step update)
c    # 2. clawpack46_time_sync_accumulate_waves (after step update)
c    # 3. clawpack46_fort_time_sync_f2c ()
c    # 4. clawpack46_fort_time_sync_copy
c    # ----------------------------------------------------------------


c    # -----------------------------------------------------------------
c    # Called BEFORE step update.  This stores the value of the flux
c    # function at cell centers.  At k=0, the flux at the the first interior
c    # interior cell is stored;  At k=1, the flux evaluated in the 
c    # ghost cell is stored;    If there is no flux function, we should 
c    # set up a dummy function that returns zero.
c    # 
c    # These will be stored for each grid and used to compute
c    # corrections.
c    # -----------------------------------------------------------------
      subroutine clawpack46_time_sync_store_flux(mx,my,mbc,meqn,
     &      maux, blockno, patchno, dt,el0, el1, el2, el3,q, aux,
     &      flux0,flux1,flux2,flux3,
     &      rpn2_cons,qvec,auxvec_center,auxvec_edge, flux)

      implicit none

      integer mx,my,mbc,meqn, maux, blockno, patchno
      double precision dt

      double precision el0(my), el1(my), el2(mx), el3(mx)

      double precision flux0(my,meqn,0:1)
      double precision flux1(my,meqn,0:1)
      double precision flux2(mx,meqn,0:1)
      double precision flux3(mx,meqn,0:1)

      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision qvec(meqn),flux(meqn)
      double precision auxvec_center(maux), auxvec_edge(maux)

      integer i,j,k,m, idir, iface

      do j = 1,my
c        # left side k=0 = in; k=1 = out
         idir = 0
         do k = 0,1
            iface = k
            do m = 1,maux
c              # Cell centered values;  Each cell-centered value
c              # will be projected onto edge between interior and
c              # exterior cell.                
               auxvec_center(m) = aux(1-k,j,m)

c              # Use this array to get edge values between interior
c              # and exterior cells.
               auxvec_edge(m) = aux(1,j,m)
            enddo


            do m = 1,meqn
               qvec(m) = q(1-k,j,m)
            enddo
            call rpn2_cons(meqn,maux,idir,iface,qvec,
     &               auxvec_center,auxvec_edge,flux)         
            do m = 1,meqn
               flux0(j,m,k) = flux0(j,m,k) + dt*el0(j)*flux(m)
            enddo
         enddo

c        # right side 0 = in; 1 = out
         do k = 0,1
            iface = 1-k
            do m = 1,maux
c              # Cell centered values                
               auxvec_center(m) = aux(mx+k,j,m)

c              # Edge between ghost cell and interior cell               
               auxvec_edge(m) = aux(mx+1,j,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(mx+k,j,m)
            enddo
            call rpn2_cons(meqn,maux,idir,iface,qvec,
     &              auxvec_center,auxvec_edge,flux)         
            do m = 1,meqn
               flux1(j,m,k) = flux1(j,m,k) + dt*el1(j)*flux(m)
            enddo
         enddo
      enddo


      idir = 1
      do i = 1,mx
c        # bottom side 0 = in; 1 = out0
         do k = 0,1
            iface = k + 2
            do m = 1,maux
c              # Cell centered values                
               auxvec_center(m) = aux(i,1-k,m)

c              # Edge between ghost cell and interior cell               
               auxvec_edge(m) = aux(i,1,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(i,1-k,m)
            enddo
            call rpn2_cons(meqn,maux,idir,iface,qvec,
     &               auxvec_center, auxvec_edge,flux)         
            do m = 1,meqn
               flux2(i,m,k) = flux2(i,m,k) + dt*el2(i)*flux(m)
            enddo
         enddo

c        # Top side 0 = in; 1 = out
         do k = 0,1
            iface = 3-k
            do m = 1,maux
c              # Cell centered values                
               auxvec_center(m) = aux(i,my+k,m)

c              # Edge between ghost cell and interior cell               
               auxvec_edge(m) = aux(i,my+1,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(i,my+k,m)
            enddo
            call rpn2_cons(meqn,maux,idir,iface,qvec,
     &              auxvec_center,auxvec_edge, flux)         
            do m = 1,meqn
               flux3(i,m,k) = flux3(i,m,k) + dt*el3(i)*flux(m)
            enddo
         enddo
      enddo

      end

c    # -----------------------------------------------------------------
c    # This is called AFTER the step update.   This accumulates plus and 
c    # minus waves, scaled by dt*edge_length.  This scaling takes care 
c    # of any division by 2 or 4. 
c    # -----------------------------------------------------------------
      subroutine clawpack46_time_sync_accumulate_waves(mx,my,
     &                                                 mbc,meqn,
     &      dt, dx,dy, patchno,
     &      el0, el1, el2, el3,
     &      fp,fm,gp,gm,
     &      fp_left,fp_right,
     &      fm_left,fm_right,
     &      gp_bottom,gp_top,
     &      gm_bottom,gm_top)

      implicit none

      integer mx,my,mbc,meqn,patchno
      double precision dt,dx,dy

      double precision fp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision fm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision fp_left(my,meqn),fp_right(my,meqn)
      double precision fm_left(my,meqn),fm_right(my,meqn)

      double precision gp_bottom(mx,meqn),gp_top(mx,meqn)
      double precision gm_bottom(mx,meqn),gm_top(mx,meqn)

      double precision el0(my), el1(my), el2(mx), el3(mx)

      integer i,j,m

c     # In flux2, fp = -apdq.  Then the update is 
c     # written as -(fm-fp) = (fp-fm) = -(apdq+amdq)
c     # We change the sign back here so that we can write
c     # the update as fp=apdq, fm = amdq.      

c     # NOTE : Fluxes have already been scaled by gamma, e.g. 
c     # x-face-length/dy or y-face-length/dx. To get scaling 
c     # by physical length, we just need to multiply by 
c     # dx or dy.


      do j = 1,my
         do m = 1,meqn
            fp_left(j,m)  = fp_left(j,m) - dt*dy*fp(1,j,m)
            fm_left(j,m)  = fm_left(j,m) + dt*dy*fm(1,j,m)

            fp_right(j,m) = fp_right(j,m)-dt*dy*fp(mx+1,j,m)
            fm_right(j,m) = fm_right(j,m)+dt*dy*fm(mx+1,j,m)
         enddo
      enddo

      do i = 1,mx
         do m = 1,meqn
            gp_bottom(i,m) = gp_bottom(i,m) - dt*dx*gp(i,1,m)
            gm_bottom(i,m) = gm_bottom(i,m) + dt*dx*gm(i,1,m)

            gp_top(i,m) = gp_top(i,m) - dt*dx*gp(i,my+1,m)
            gm_top(i,m) = gm_top(i,m) + dt*dx*gm(i,my+1,m)
         enddo
      enddo

      return

      end

c    # -----------------------------------------------------------------
c    # Add fine grid corrections to coarse grid.  
c    # -----------------------------------------------------------------

      subroutine clawpack46_fort_time_sync_f2c (mx,my,mbc,meqn,
     &                                          idir, iface_coarse,
     &                                          coarse_blockno,
     &                                          fine_blockno,
     &                                          normal_match,
     &                                          area0, area1,
     &                                          area2, area3,
     &                                          qcoarse,
     &                                          fpthis0,fmthis1,
     &                                          gpthis2,gmthis3,
     &                                          fmfine0, fpfine1,
     &                                          gmfine2, gpfine3,
     &                                          efthis0, efthis1, 
     &                                          efthis2, efthis3,
     &                                          eff0, eff1, eff2, eff3,
     &                                          qfine_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse
      integer coarse_blockno, fine_blockno, normal_match

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # Double number of ghost cells so they map to whole coarse
c     # grid cells      
      double precision qfine_dummy(1-2*mbc:mx+2*mbc,
     &                              1-2*mbc:my+2*mbc,meqn)

      double precision area0(my), area1(my), area2(mx), area3(mx)

      double precision fpthis0(my,meqn)
      double precision fmthis1(my,meqn)
      double precision gpthis2(mx,meqn)
      double precision gmthis3(mx,meqn)

      double precision fmfine0(my,meqn)
      double precision fpfine1(my,meqn)
      double precision gmfine2(mx,meqn)
      double precision gpfine3(mx,meqn)

      double precision eff0(my,meqn,0:1)
      double precision eff1(my,meqn,0:1)
      double precision eff2(mx,meqn,0:1)
      double precision eff3(mx,meqn,0:1)

      double precision efthis0(my,meqn,0:1)
      double precision efthis1(my,meqn,0:1)
      double precision efthis2(mx,meqn,0:1)
      double precision efthis3(mx,meqn,0:1)

      integer*8 transform_cptr

      logical is_valid_correct

      double precision fm,fp,gm,gp,ef

      integer i, j, ic, jc, ii1, ii2, jj1, jj2, ii, jj
      integer i2(0:3), j2(0:3), m, mq
      double precision deltac, areac


      integer a(2,2), f(2), sc, nm

      logical is_valid_average, skip_this_grid

      call fclaw2d_clawpatch_build_transform(transform_cptr,a,f)

c     # This is here to test the normal_match value; eventually, it
c     # should go away.      

      if (idir .eq. 0) then
        sc = (a(1,1) + a(2,1))/2
      else
        sc = (a(1,2) + a(2,2))/2
      endif

c     # Ignore call to normal_match in calling routine      
      if (normal_match .eq. 0) then
          nm = -1
      else
          nm = 1          
      endif
c      nm = sc 
      if (sc .ne. nm) then
          write(6,*) 'time_sync : Normal match does not match'
          write(6,*) 'sc = ', sc
          write(6,*) 'normal_match = ', nm
          write(6,*) a(1,1), a(2,1), a(1,2),a(2,2)
          stop
      endif


c     # We do not know which edge will be used, or how it will be used
c     # so we fill in data at every edge.  Use double the number of 
c     # ghost cells so that data aligns with coarse grid ghost cells.
      do mq = 1,meqn
         do j = 1,my/2
            jj1 = 2*(j-1) + 1
            jj2 = 2*(j-1) + 2
            fm = fmfine0(jj1,mq) + fmfine0(jj2,mq)
            ef = eff0(jj1,mq,1) + eff0(jj2,mq,1)
            do jj = jj1,jj2
                do ii = -1,0
                    qfine_dummy(ii,jj,mq)  = fm
                    qfine_dummy(ii-2,jj,mq) = ef
                enddo
            enddo
     
            fp = fpfine1(jj1,mq) + fpfine1(jj2,mq)
            ef = eff1(jj1,mq,1) + eff1(jj2,mq,1)
            do jj = jj1,jj2
                do ii = 1,2
                    qfine_dummy(mx+ii,jj,mq) = fp
                    qfine_dummy(mx+ii+2,jj,mq) = ef
                enddo
            enddo
         enddo

         do i = 1,mx/2
            ii1 = 2*(i-1) + 1
            ii2 = 2*(i-1) + 2
            gm = gmfine2(ii1,mq) + gmfine2(ii2,mq)
            ef = eff2(ii1,mq,1) + eff2(ii2,mq,1)
            do ii = ii1,ii2
                do jj = -1,0
                    qfine_dummy(ii,jj,mq)  = gm
                    qfine_dummy(ii,jj-2,mq) = ef
                enddo
            enddo
            gp = gpfine3(ii1,mq) + gpfine3(ii2,mq)
            ef = eff3(ii1,mq,1) + eff3(ii2,mq,1)
            do ii = ii1,ii2
                do jj = 1,2
                    qfine_dummy(ii,my+jj,mq) = gp
                    qfine_dummy(ii,my+jj+2,mq) = ef
                enddo
            enddo

         enddo
      enddo

      if (idir .eq. 0) then
c         # sign change ('sc') to account for normals at 0-1
c         # patch faces which may point in a different direction.
c         # First column is A*[1;0]        
          sc = (a(1,1) + a(2,1))/2
          do mq = 1,meqn
              do jc = 1,my
                  if (iface_coarse .eq. 0) then
                      ic = 1
                      call fclaw2d_clawpatch_transform_face_half(ic,jc,
     &                     i2,j2,transform_cptr)
                      deltac = 0
                      areac = area0(jc)
                      if (.not. skip_this_grid(idir,i2,j2,mx,my)) then
                          fp = qfine_dummy(i2(0),j2(0),mq)
                          call 
     &                fclaw2d_clawpatch_transform_face_half(ic+1,jc,
     &                         i2,j2,transform_cptr)
c                         # Check for validity?                          
                          ef = sc*qfine_dummy(i2(0),j2(0),mq)
                          deltac = (ef-efthis0(jc,mq,0))-
     &                            (fp-fpthis0(jc,mq))
                      endif
                  elseif (iface_coarse .eq. 1) then
                      ic = mx
                      call 
     &         fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &                      transform_cptr)
                      deltac = 0
                      areac = area1(jc)
                      if (.not. skip_this_grid(idir,i2,j2,mx,my)) then
                          fm = qfine_dummy(i2(0),j2(0),mq)
    
                          call 
     &       fclaw2d_clawpatch_transform_face_half(ic-1,jc,
     &                          i2,j2,transform_cptr)
                          ef = sc*qfine_dummy(i2(0),j2(0),mq)
    
                          deltac = (efthis1(jc,mq,0)-ef)-
     &                             (fm-fmthis1(jc,mq))
                      endif
                  endif
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) + deltac/areac
              enddo
          enddo
      else
c         # sign change ('sc') to account for normals at 2-3 
c         # patch faces which may point in a different direction.
c         # Second column is A*[0;1]        
          sc = (a(1,2) + a(2,2))/2
          do mq = 1,meqn
              do ic = 1,mx                  
                  if (iface_coarse .eq. 2) then
                      jc = 1
                      call 
     &     fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &                      transform_cptr)
                      deltac = 0
                      areac = area2(ic)
                      if (.not. skip_this_grid(idir,i2,j2,mx,my)) then                
                          gp = qfine_dummy(i2(0),j2(0),mq)
                          call 
     &      fclaw2d_clawpatch_transform_face_half(ic,jc+1,
     &                         i2,j2,transform_cptr)
                          ef = sc*qfine_dummy(i2(0),j2(0),mq)
    
                          deltac = (ef-efthis2(ic,mq,0))-
     &                        (gp-gpthis2(ic,mq))
                      endif
                  else if (iface_coarse .eq. 3) then
                      jc = my
                      call 
     &       fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2,
     &                      transform_cptr)
                      deltac = 0
                      areac = area3(ic)
                      if (.not. skip_this_grid(idir,i2,j2,mx,my)) then                   
                          gm = qfine_dummy(i2(0),j2(0),mq)
                          call 
     &         fclaw2d_clawpatch_transform_face_half(ic,jc-1,
     &                          i2,j2,transform_cptr)
                          ef = sc*qfine_dummy(i2(0),j2(0),mq)
    
                         deltac = (efthis3(ic,mq,0)-ef)-
     &                            (gm-gmthis3(ic,mq))
                      endif
                  endif
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) + deltac/areac
              enddo
          enddo
      endif

      end


c    # -----------------------------------------------------------------
c    # Add wave corrections at same level interfaces.  This accounts for
c    # metric mismatches that can occur at block boundaries.
c    # -----------------------------------------------------------------
      subroutine clawpack46_fort_time_sync_samesize(mx,my,mbc,meqn,
     &                                          idir, this_iface,
     &                                          this_blockno, 
     &                                          neighbor_blockno,
     &                                          area0, area1,
     &                                          area2, area3,
     &                                          qthis,
     &                                          fpthis0,fmthis1,
     &                                          gpthis2,gmthis3,
     &                                          fmneighbor0, 
     &                                          fpneighbor1,
     &                                          gmneighbor2, 
     &                                          gpneighbor3,
     &                                          efthis0, 
     &                                          efthis1, 
     &                                          efthis2, 
     &                                          efthis3,
     &                                          efneighbor0, 
     &                                          efneighbor1, 
     &                                          efneighbor2, 
     &                                          efneighbor3,
     &                                          qthis_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse
      integer this_iface, this_blockno, neighbor_blockno

      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision qthis_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision area0(my), area1(my), area2(mx), area3(mx)

      double precision fpthis0(my,meqn)
      double precision fmthis1(my,meqn)
      double precision gpthis2(mx,meqn)
      double precision gmthis3(mx,meqn)

      double precision fmneighbor0(my,meqn)
      double precision fpneighbor1(my,meqn)
      double precision gmneighbor2(mx,meqn)
      double precision gpneighbor3(mx,meqn)

c     # These store flux evaluations f(q) at the first layer of ghost 
c     # cells and the first interior row/column.  Interior values are 
c     # stored at the '0' index; ghost layer stored at the '1' index.    
      double precision efthis0(my,meqn,0:1)
      double precision efthis1(my,meqn,0:1)
      double precision efthis2(mx,meqn,0:1)
      double precision efthis3(mx,meqn,0:1)

      double precision efneighbor0(my,meqn,0:1)
      double precision efneighbor1(my,meqn,0:1)
      double precision efneighbor2(mx,meqn,0:1)
      double precision efneighbor3(mx,meqn,0:1)


      double precision delta, deltap, deltam, area

      integer a(2,2), f(2), sc

      double precision fp,fm,gp,gm,ef

      double precision neighborval

      integer*8 transform_cptr


      integer i,j,mq
      integer i1,j1, i2, j2


      call build_transform_samesize(transform_cptr,a,f)

      idir = this_iface/2

      do mq = 1,meqn
         do j = 1,my
             qthis_dummy(0,j,mq)  = fmneighbor0(j,mq)
             qthis_dummy(-1,j,mq) = efneighbor0(j,mq,1)

             qthis_dummy(mx+1,j,mq) = fpneighbor1(j,mq)
             qthis_dummy(mx+2,j,mq) = efneighbor1(j,mq,1)
         enddo

         do i = 1,mx
             qthis_dummy(i,0,mq)  = gmneighbor2(i,mq)
             qthis_dummy(i,-1,mq) = efneighbor2(i,mq,1)

             qthis_dummy(i,my+1,mq) = gpneighbor3(i,mq)
             qthis_dummy(i,my+2,mq) = efneighbor3(i,mq,1)
         enddo
      enddo

      if (idir .eq. 0) then
c         # Get sign change to account for possible mismatch in normals    
          sc = a(1,1) + a(2,1)    
          do mq = 1,meqn
              do j = 1,my
c                 # Exchange at low side of 'this' grid in
c                 # x-direction (idir == 0) 
                  j1 = j
                  if (this_iface .eq. 0) then
                      i1 = 1
                      call 
     &       fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,
     &                      transform_cptr)
                      fp = qthis_dummy(i2,j2,mq)
                      call 
     &       fclaw2d_clawpatch_transform_face(i1+1,j1,i2,j2,
     &                      transform_cptr)
                      ef = sc*qthis_dummy(i2,j2,mq)

                      delta = (ef-efthis0(j1,mq,0))-(fp-fpthis0(j1,mq))
                      area = area0(j1)

                  else
                      i1 = mx
                      call 
     &       fclaw2d_clawpatch_transform_face(i1,j1,i2,j2,
     &                      transform_cptr)
                      fm = qthis_dummy(i2,j2,mq)
                      call 
     &       fclaw2d_clawpatch_transform_face(i1-1,j1,i2,j2,
     &                      transform_cptr)
                      ef = sc*qthis_dummy(i2,j2,mq)

                      delta = (efthis1(j1,mq,0)-ef)-(fm-fmthis1(j1,mq))
                      area = area1(j1)

                  endif
                  qthis(i1,j1,mq) = qthis(i1,j1,mq) + 0.5*delta/area
              enddo
          enddo
      else
c         # Get sign change ('sc') to account for possible 
c         # mismatch in direction of normal vectors
          sc = a(1,2) + a(2,2)    
          do mq = 1,meqn
              do i = 1,mx
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  i1 = i
                  if (this_iface .eq. 2) then
                      j1 = 1
                      call fclaw2d_clawpatch_transform_face(i1,j1,
     &                 i2,j2,transform_cptr)
                      gp = qthis_dummy(i2,j2,mq)
                      call fclaw2d_clawpatch_transform_face(i1,j1+1,
     &                      i2,j2,transform_cptr)
                      ef = sc*qthis_dummy(i2,j2,mq)

                      delta = (ef-efthis2(i1,mq,0))-(gp-gpthis2(i1,mq))
                      area = area2(i1)

                  else
c                     # Coarse grid is bottom grid; fine is top grid   
c                     # efc3(0) is flux stored in the coarse grid 
c                     # interior cell.    
                      j1 = my
                      call fclaw2d_clawpatch_transform_face(i1,j1,
     &                        i2,j2,transform_cptr)
c                      write(6,*) this_iface, i1, j1, i2, j2
                      gm = qthis_dummy(i2,j2,mq)
                      call fclaw2d_clawpatch_transform_face(i1,j1-1,
     &                         i2,j2,transform_cptr)
                      ef = sc*qthis_dummy(i2,j2,mq)

                     delta = (efthis3(i1,mq,0)-ef)-(gm-gmthis3(i1,mq))
                     area = area3(i1)

                  endif
                  qthis(i1,j1,mq) = qthis(i1,j1,mq) + 0.5*delta/area                 
              enddo
          enddo
      endif
      end

      logical function skip_this_grid(idir,i2,j2,mx,my)
      implicit none

      integer idir, i2(0:3), j2(0:3), mx,my
      integer m
      logical is_valid_correct

      skip_this_grid = .false.
      do m = 0,3
          if (.not. is_valid_correct(idir,i2(m),j2(m),mx,my)) then
              skip_this_grid = .true.
              exit
          endif
      enddo   


      end


      logical function is_valid_correct(idir,i,j,mx,my)
      implicit none

      integer i,j,mx,my, idir
      logical i1, j1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my

c     # At least one of the above should be true.
      is_valid_correct = xor(i1,j1)


      end


      subroutine build_transform_samesize(transform_ptr,a,f)
      implicit none

      integer a(2,2)
      integer*8 transform_ptr
      integer f(2)
      integer icn, jcn
      integer i1,j1

c     # Assume index mapping fclaw2d_transform_face has the
c     # the form
c     #
c     #       T(ic,jc) = A*(ic,jc) + F = (icn,jcn)
c     #
c     # where (ic,jc) is the coarse grid index, and (icn,jcn)
c     # is the neighbor grid index.
c     #
c     # We can recover coefficients A(2,2) with the following
c     # calls to T.

      i1 = 0
      j1 = 0
      call fclaw2d_clawpatch_transform_face(i1,j1,icn,jcn,
     &      transform_ptr)
      f(1) = icn
      f(2) = jcn

      i1 = 1
      j1 = 0
      call fclaw2d_clawpatch_transform_face(i1,j1,icn,jcn,
     &      transform_ptr)
      a(1,1) = icn - f(1)
      a(2,1) = jcn - f(2)

      i1 = 0
      j1 = 1
      call fclaw2d_clawpatch_transform_face(i1,j1,icn,jcn,
     &      transform_ptr)
      a(1,2) = icn - f(1)
      a(2,2) = jcn - f(2)

      end


