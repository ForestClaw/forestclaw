c    # ----------------------------------------------------------------
c    # This file contains routines which accumulate fluxes, waves 
c    # and add corrections to coarse grid at edges of both same size
c    # grids and coarse/fine interfaces. 
c    # 
c    # 1. clawpack46_cons_update_store_flux (before step update)
c    # 2. clawpack46_cons_update_accumulate_waves (after step update)
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
      subroutine clawpack46_update_cons_store_flux(mx,my,mbc,meqn,
     &      maux, dt,el0, el1, el2, el3,q, aux,
     &      flux0,flux1,flux2,flux3,
     &      rpn2_cons,qvec,auxvec,flux)

      implicit none

      integer mx,my,mbc,meqn, maux,idir
      double precision dt

      double precision el0(my), el1(my), el2(mx), el3(mx)

      double precision flux0(my,meqn,0:1)
      double precision flux1(my,meqn,0:1)
      double precision flux2(mx,meqn,0:1)
      double precision flux3(mx,meqn,0:1)

      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision qvec(meqn),auxvec(maux),flux(meqn)

      integer i,j,k,m

      do j = 1,my
c        # left side k=0 = in; k=1 = out
        idir = 0
         do k = 0,1
            do m = 1,maux
               auxvec(m) = aux(1-k,j,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(1-k,j,m)
            enddo
            call rpn2_cons(meqn,maux,idir,qvec,auxvec,flux)         
            do m = 1,meqn
               flux0(j,m,k) = flux0(j,m,k) + dt*el0(j)*flux(m)
            enddo
         enddo

c        # right side 0 = in; 1 = out
         do k = 0,1
            do m = 1,maux
               auxvec(m) = aux(mx+k,j,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(mx+k,j,m)
            enddo
            call rpn2_cons(meqn,maux,idir,qvec,auxvec,flux)         
            do m = 1,meqn
               flux1(j,m,k) = flux1(j,m,k) + dt*el1(j)*flux(m)
            enddo
         enddo
      enddo


      do i = 1,mx
c        # bottom side 0 = in; 1 = out
         idir = 1
         do k = 0,1
            do m = 1,maux
               auxvec(m) = aux(i,1-k,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(i,1-k,m)
            enddo
            call rpn2_cons(meqn,maux,idir,qvec,auxvec,flux)         
            do m = 1,meqn
               flux2(i,m,k) = flux2(i,m,k) + dt*el2(i)*flux(m)
            enddo
         enddo

c        # Top side 0 = in; 1 = out
         do k = 0,1
            do m = 1,maux
               auxvec(m) = aux(i,my+k,m)
            enddo
            do m = 1,meqn
               qvec(m) = q(i,my+k,m)
            enddo
            call rpn2_cons(meqn,maux,idir,qvec,auxvec,flux)         
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
      subroutine clawpack46_cons_update_accumulate_waves(mx,my,
     &                                                   mbc,meqn,
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

      do j = 1,my
         do m = 1,meqn
c           # In flux2, fp = -apdq.  Then the update is 
c           # written as -(fm-fp) = (fp-fm) = -(apdq+amdq)
c           # We change the sign back here so that we can write
c           # the update as fp=apdq, fm = amdq.      

c           # NOTE : Fluxes have already been scaled by gamma, e.g. 
c           # x-face-length/dy or y-face-length/dx

            fp_left(j,m) = fp_left(j,m) - dt*dy*fp(1,j,m)
            fm_left(j,m) = fm_left(j,m) + dt*dy*fm(1,j,m)

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
     &                                          area0, area1,
     &                                          area2, area3,
     &                                          qcoarse,
     &                                          fpcoarse0,fmcoarse1,
     &                                          gpcoarse2,gmcoarse3,
     &                                          fmfine0, fpfine1,
     &                                          gmfine2, gpfine3,
     &                                          efc0, efc1, efc2, efc3,
     &                                          eff0, eff1, eff2, eff3,
     &                                          maskfine, qfine_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision qfine_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer maskfine(1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision area0(my), area1(my), area2(mx), area3(mx)

      double precision fpcoarse0(my,meqn)
      double precision fmcoarse1(my,meqn)
      double precision gpcoarse2(mx,meqn)
      double precision gmcoarse3(mx,meqn)

      double precision fmfine0(my,meqn)
      double precision fpfine1(my,meqn)
      double precision gmfine2(mx,meqn)
      double precision gpfine3(mx,meqn)

      double precision eff0(my,meqn,0:1)
      double precision eff1(my,meqn,0:1)
      double precision eff2(mx,meqn,0:1)
      double precision eff3(mx,meqn,0:1)

      double precision efc0(my,meqn,0:1)
      double precision efc1(my,meqn,0:1)
      double precision efc2(mx,meqn,0:1)
      double precision efc3(mx,meqn,0:1)

      integer*8 transform_cptr

      integer i,j,ii1,ii2,jj1,jj2,ii,jj, m, mq
      integer ic,jc, r2, refratio
      double precision deltam,deltap, deltac, sum,areac, fineval
      logical is_valid_correct

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      logical is_valid_average, skip_this_grid


      refratio = 2
      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'cons_coarse_correct ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


c    # Visit each fine grid face and fill in correction term into
c    # qfine_dummy. We don't really know which face corresponds
c    # to the common face shared between fine and coarse grids, so we
c    # do them all.  Then use the transform to select correct the 
c    # fine grid correction for the coarse grid.

      do mq = 1,meqn
         do j = 1,my/2
            jj1 = 2*(j-1) + 1     !! indexing starts at 1
            jj2 = 2*(j-1) + 2
            deltam = fmfine0(jj1,mq) + fmfine0(jj2,mq) +
     &               eff0(jj1,mq,1) + eff0(jj2,mq,1)
     
            deltap = (fpfine1(jj1,mq) + fpfine1(jj2,mq)) -
     &               (eff1(jj1,mq,1) + eff1(jj2,mq,1))

            if (j .eq. 1) then
c                write(6,200) fmfine0(jj1,mq),eff0(jj1,mq,1)
            endif
 200        format(2F16.8)

c           # Put the same value in four ghost cells;  grab the first
c           # one that shows up in the transformation.          
            do jj = jj1,jj2
               do ii = 1,mbc
                  qfine_dummy(1-ii,jj,mq) = -deltam
                  qfine_dummy(mx+ii,jj,mq) =-deltap
               enddo
            enddo
         enddo

         do i = 1,mx/2
            ii1 = 2*(i-1) + 1        !! indexing starts at 1
            ii2 = 2*(i-1) + 2
            deltam = gmfine2(ii1,mq) + gmfine2(ii2,mq) +
     &               eff2(ii1,mq,1) + eff2(ii2,mq,1)

            deltap = gpfine3(ii1,mq) + gpfine3(ii2,mq) -
     &               (eff3(ii1,mq,1) + eff3(ii2,mq,1))

            do ii = ii1,ii2
               do jj = 1,mbc
                  qfine_dummy(ii,1-jj,mq) = -deltam
                  qfine_dummy(ii,my+jj,mq) = -deltap
               enddo
            enddo
         enddo
      enddo


c     # Add corrections from above to coarse grid.  Use transforms to get
c     # correct location.
      if (idir .eq. 0) then
         do mq = 1,meqn
            do jc = 1,my
              
c              # Move this to beginning of routine.              
               if (iface_coarse .eq. 0) then
                  ic = 1
               elseif (iface_coarse .eq. 1) then
                  ic = mx
               endif

               call fclaw2d_transform_face_half(ic,jc,i2,j2,
     &               transform_cptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. is_valid_correct(idir,i2(m),j2(m),mx,my))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo

               if (.not. skip_this_grid) then
                  fineval = qfine_dummy(i2(1),j2(1),mq)
                  if (iface_coarse .eq. 0) then
c                    # Coarse grid is right; fine grid is left
c                    # efc0(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + fpcoarse0(jc,mq) - efc0(jc,mq,0)
                     areac = area0(jc)

c                    # Reset these, since they will no longer be needed                     
                     fpcoarse0(jc,mq) = 0
                     efc0(jc,mq,0) = 0
                     efc0(jc,mq,1) = 0

                  else
c                   # Coarse grid is left; fine grid is right                    
c                    # efc1(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + fmcoarse1(jc,mq) + efc1(jc,mq,0)
                     areac = area1(jc)

c                    # Reset these, since they will no longer be needed                     
                     fmcoarse1(jc,mq) = 0
                     efc1(jc,mq,0) = 0
                     efc1(jc,mq,1) = 0
                  endif   
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) + deltac/areac
               endif
            enddo
         enddo
      else
         do mq = 1,meqn
            do ic = 1,mx
               if (iface_coarse .eq. 2) then
                  jc = 1
               elseif (iface_coarse .eq. 3) then
                  jc = my
               endif

               call fclaw2d_transform_face_half(ic,jc,i2,j2,
     &               transform_cptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. is_valid_correct(idir,i2(m),j2(m),mx,my))
     &                  then
                     skip_this_grid = .true.
                  endif
               enddo
               if (.not. skip_this_grid) then
                  fineval = qfine_dummy(i2(1),j2(1),mq)
                  if (iface_coarse .eq. 2) then
c                    # Coarse is top grid  Fine grid is bottom grid; 
c                    # efc2(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + gpcoarse2(ic,mq) - efc2(ic,mq,0)
                     areac = area2(ic)

c                    # Reset these, since they will no longer be needed                     
                     gpcoarse2(ic,mq) = 0
                     efc2(ic,mq,0) = 0
                     efc2(ic,mq,1) = 0

                  else
c                    # Coarse grid is bottom grid; fine is top grid   
c                    # efc3(0) is flux stored in the coarse grid 
c                    # interior cell.                   
                     deltac = fineval + gmcoarse3(ic,mq) + efc3(ic,mq,0)
                     areac = area3(ic)

c                    # Reset these, since they will no longer be needed                     
                     gmcoarse3(ic,mq) = 0
                     efc3(ic,mq,0) = 0
                     efc3(ic,mq,1) = 0

                  endif
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) + deltac/areac
               endif                    !! skip grid loop
            enddo
         enddo
      endif


      end


c    # -----------------------------------------------------------------
c    # Add wave corrections at same level interfaces.  This accounts for
c    # metric mismatches that can occur at block boundaries.
c    # -----------------------------------------------------------------
      subroutine clawpack46_fort_time_sync_copy(mx,my,mbc,meqn,
     &                                          idir, this_iface,
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
     &                                          mask, qneighbor_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse
      integer this_iface

      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision qneighbor_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer mask(1-mbc:mx+mbc,1-mbc:my+mbc)

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


      double precision delta, deltap, deltam, areathis

      double precision neighborval

      integer*8 transform_cptr


      integer i,j,ibc,jbc,mq
      integer i1,j1, i2, j2

      idir = this_iface/2


      do mq = 1,meqn
         do j = 1,my
c           # Left edge - "this" grid is right grid; 
c           # "neighbor" grid is left grid.          
            deltam = fmneighbor0(j,mq)  + efneighbor0(j,mq,1)

c           # Right edge                         
            deltap = fpneighbor1(j,mq) - efneighbor1(j,mq,1)

            qneighbor_dummy(0,j,mq) = -deltam
            qneighbor_dummy(mx+1,j,mq) = -deltap
         enddo

         do i = 1,mx
            deltam = gmneighbor2(i,mq) + efneighbor2(i,mq,1)
            deltap = gpneighbor3(i,mq) - efneighbor3(i,mq,1)

            qneighbor_dummy(i,0,mq)    = -deltam
            qneighbor_dummy(i,my+1,mq) = -deltap
         enddo
      enddo



c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      do mq = 1,meqn
         if (idir .eq. 0) then
            do j = 1,my
               do ibc = 1,mbc
c                 # Exchange at low side of 'this' grid in
c                 # x-direction (idir == 0) 
                  j1 = j
                  if (this_iface .eq. 0) then
                     i1 = 1
                  elseif (this_iface .eq. 1) then
                     i1 = mx
                  endif
                  call fclaw2d_transform_face(i1,j1,i2,j2,
     &                  transform_cptr)

                  neighborval = qneighbor_dummy(i2,j2,mq)
                  if (this_iface .eq. 0) then
                     delta = neighborval  + fpthis0(j1,mq) 
     &                             - efthis0(j1,mq,0)
                     areathis = area0(j1)

                     fpthis0(j1,mq) = 0
                     efthis0(j1,mq,0) = 0
                     efthis0(j1,mq,1) = 0
                  else
c                   # Coarse grid is left; fine grid is right                    
c                    # efc1(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     delta = neighborval + fmthis1(j1,mq) + 
     &                             efthis1(j1,mq,0)
                     areathis = area1(j1)

c                    # Reset these, since they will no longer be needed                     
                     fmthis1(j1,mq) = 0
                     efthis1(j1,mq,0) = 0
                     efthis1(j1,mq,1) = 0
                  endif
                  qthis(i1,j1,mq) = qthis(i1,j1,mq) + delta/areathis
               enddo
            enddo
         else
            do jbc = 1,mbc
               do i = 1,mx
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  i1 = i
                  if (this_iface .eq. 2) then
                     j1 = 1
                  elseif (this_iface .eq. 3) then
                     j1 = my
                  endif
                  call fclaw2d_transform_face(i1,j1,i2,j2,
     &                  transform_cptr)
                  neighborval = qneighbor_dummy(i2,j2,mq)
                  if (this_iface .eq. 2) then
c                    # Coarse is top grid  Fine grid is bottom grid; 
c                    # efc2(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     delta = neighborval + gpthis2(i1,mq) 
     &                          - efthis2(i1,mq,0)
                     areathis = area2(i1)

c                    # Reset these, since they will no longer be needed                     
                     gpthis2(i1,mq) = 0
                     efthis2(i1,mq,0) = 0
                     efthis2(i1,mq,1) = 0
                  else
c                    # Coarse grid is bottom grid; fine is top grid   
c                    # efc3(0) is flux stored in the coarse grid 
c                    # interior cell.                   
                     delta = neighborval + gmthis3(i1,mq) + 
     &                            efthis3(i1,mq,0)
                     areathis = area3(i1)

c                    # Reset these, since they will no longer be needed                     
                     gmthis3(i1,mq) = 0
                     efthis3(i1,mq,0) = 0
                     efthis3(i1,mq,1) = 0
                  endif
                  qthis(i1,j1,mq) = qthis(i1,j1,mq) + delta/areathis                 
               enddo
            enddo
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

