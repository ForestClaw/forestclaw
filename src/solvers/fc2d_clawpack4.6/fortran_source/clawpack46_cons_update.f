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
c    # function at cell centers.  At k=0, the flux evaluted in the 
c    # ghost cell is stored;  at k=1, the flux at the the first interior
c    # cell is stored;  If there is no flux function, we should set up a 
c    # dummy function that returns zero.
c    # -----------------------------------------------------------------
      subroutine clawpack46_cons_update_store_flux(mx,my,mbc,meqn,
     &      maux, dt,
     &      el0, el1, el2, el3,
     &      q, aux,
     &      flux0,flux1,flux2,flux3,
     &      rpn2_cons,qvec,auxvec,flux)

      implicit none

      integer mx,my,mbc,meqn, maux,idir
      double precision dt

      double precision el0(my), el1(my), el2(mx), el3(mx)

      double precision flux0(my,meqn,0:1)
      double precision flux1(my,meqn,0:1)
      double precision flux2(my,meqn,0:1)
      double precision flux3(my,meqn,0:1)

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
               flux0(j,m,k) = dt*el0(j)*flux(m)
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
               flux1(j,m,k) = dt*el1(j)*flux(m)
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
               flux2(i,m,k) = dt*el2(i)*flux(m)
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
               flux3(i,m,k) = dt*el3(i)*flux(m)
            enddo
         enddo
      enddo

      end

c    # -----------------------------------------------------------------
c    # This is called AFTER the step update.   This accumulates plus and 
c    # minus waves, scaled by dt*edge_length.  
c    # -----------------------------------------------------------------
      subroutine clawpack46_cons_update_accumulate_wave(mx,my,
                                                        mbc,meqn,
     &      dt, patchno,
     &      el0, el1, el2, el3,
     &      fp,fm,gp,gm,
     &      fp_left,fp_right,
     &      fm_left,fm_right,
     &      gp_bottom,gp_top,
     &      gm_bottom,gm_top)

      implicit none

      integer mx,my,mbc,meqn,patchno
      double precision dt

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
            fp_left(j,m) = fp_left(j,m) - dt*el0(j)*fp(1,j,m)
            fm_left(j,m) = fm_left(j,m) + dt*el0(j)*fm(1,j,m)

            fp_right(j,m) = fp_right(j,m) - dt*el1(j)*fp(mx+1,j,m)
            fm_right(j,m) = fm_right(j,m) + dt*el1(j)*fm(mx+1,j,m)
         enddo
      enddo

      do i = 1,mx
         do m = 1,meqn
            gp_bottom(i,m) = gp_bottom(i,m) - dt*el2(i)*gp(i,1,m)
            gm_bottom(i,m) = gm_bottom(i,m) + dt*el2(i)*gm(i,1,m)

            gp_top(i,m) = gp_top(i,m) - dt*el3(i)*gp(i,my+1,m)
            gm_top(i,m) = gm_top(i,m) + dt*el3(i)*gm(i,my+1,m)
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
            deltap = fpfine1(jj1,mq) + fpfine1(jj2,mq) +
     &               eff1(jj1,mq,1) + eff1(jj2,mq,1)

            do jj = jj1,jj2
               do ii = 1,mbc
                  qfine_dummy(1-ii,jj,mq) = -deltam
                  qfine_dummy(mx+ii,jj,mq) = -deltap
               enddo
            enddo
         enddo

         do i = 1,mx/2
            ii1 = 2*(i-1) + 1        !! indexing starts at 1
            ii2 = 2*(i-1) + 2
            deltam = gmfine2(ii1,mq) + gmfine3(ii2,mq) +
     &               eff2(ii1,mq,1) + eff2(ii2,mq,1)
            deltap = gpfine1(ii1,mq) + gpfine1(ii2,mq) +
     &               eff3(ii1,mq,1) + eff3(ii2,mq,1)
            do ii = ii1,ii2
               do jj = 1,mbc
                  qfine_dummy(ii,1-jj,mq) = -deltam
                  qfine_dummy(ii,my+jj,mq) = -deltap
               enddo
            enddo
         enddo
      enddo


c     # Average fine grid onto coarse grid
      if (idir .eq. 0) then
         do mq = 1,meqn
            do jc = 1,my
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
c                    # Coarse face is right; fine grid is left
c                    # efc0(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + fpcoarse0(jc,mq) + efc0(jc,mq,0)
                     areac = area0(jc)
                  else
c                    # Coarse face is left; fine grid is right                    
c                    # efc1(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + fmcoarse1(jc,mq) + efc1(jc,mq,0)
                     areac = area1(jc)
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
c                    # Fine grid is bottom grid; coarse is top grid  
c                    # efc2(0) is flux stored in the coarse grid 
c                    # interior cell. 
                     deltac = fineval + gpcoarse2(ic,mq) + efc2(ic,mq,0)
                     areac = area2(ic)
                  else
c                    # coarse grid is bottom grid; fine is top grid   
c                    # efc3(0) is flux stored in the coarse grid 
c                    # interior cell.                   
                     deltac = fineval + gmcoarse3(ic,mq) + efc3(ic,mq,0)
                     areac = area3(ic)
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
     &                                          fmneighbor0, fpneighbor1,
     &                                          gmneighbor2, gpneighbor3,
     &                                          efc0, efc1, efc2, efc3,
     &                                          eff0, eff1, eff2, eff3,
     &                                          mask, qneighbor_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse

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


      integer*8 transform_cptr


      integer i,j,ibc,jbc,mq, idir
      integer i1,j1, i2, j2

      idir = iface/2


      do mq = 1,meqn
         do j = 1,my
c           # Left edge - "this" grid is right grid; 
c           # "neighbor" grid is left grid.          
            deltam = fmneighbor0(j,mq)  - efneighbor0(j,mq,1)

c           # Right edge                         
            deltap = fpneighbor1(j,mq) - efneighbor1(j,mq,1)

            qneighbor_dummy(0,j,mq) = -deltam
            qneighbor_dummy(mx+1,j,mq) = -deltap
         enddo

         do i = 1,mx
            deltam = gmneighbor2(i,mq) + efneighbor2(i,mq,1)
            deltap = gpneighbor3(i,mq) + efneighbor3(i,mq,1)

            qneighbor_dummy(i,0,mq) = -deltam
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
                  if (iface .eq. 0) then
                     i1 = 1-ibc
                     j1 = j
                  elseif (iface .eq. 1) then
                     i1 = mx+ibc
                     j1 = j
                  endif
                  call fclaw2d_transform_face(i1,j1,i2,j2,
     &                  transform_ptr)

                  neighborval = qneighbor_dummy(i2,j2,mq)
                  if (iface_this .eq. 0) then
                     delta = neighborval + fm0()
                     qthis(i1,j1,mq) = qthis(i1,j1,mq) + delta
                  endif

               enddo
            enddo
         else
            do jbc = 1,mbc
               do i = 1,mx
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  if (iface .eq. 2) then
                     i1 = i
                     j1 = 1-jbc
                  elseif (iface .eq. 3) then
                     i1 = i
                     j1 = my+jbc
                  endif
                  call fclaw2d_transform_face(i1,j1,i2,j2,
     &                  transform_ptr)
                  qthis(i1,j1,mq) = qneighbor(i2,j2,mq)

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

      is_valid_correct = xor(i1,j1)


      end

c    # --------------------------------------------------------
c    # THIS IS NOT USED .... AND SHOULD BE DELETED!
c    # --------------------------------------------------------
c    # Called from interpolation routine (with possible context
c    # switching). This fills in coarse grid values qc/auxc
c    # on fine grids.
      subroutine clawpack46_fort_cons_coarse_to_fine(mx,my,mbc,maux,
     &       meqn, maskfine, qcoarse,auxcoarse,
     &       qfine_dummy,auxfine_dummy,
     &       idir,iface_coarse,qc0,qc1,qc2,qc3,
     &       auxc0,auxc1,auxc2, auxc3,transform_ptr)

      implicit none

      integer mx,my,mbc,maux,meqn,idir,iface_coarse
      integer maskfine(1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision auxfine_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision qc0(my,meqn), auxc0(my,maux)
      double precision qc1(my,meqn), auxc1(my,maux)
      double precision qc2(mx,meqn), auxc2(mx,maux)
      double precision qc3(mx,meqn), auxc3(mx,maux)

      integer*8 transform_ptr

c     #  Hard-coding refratio = 2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      logical is_valid_correct
      logical skip_this_grid

      integer a(2,2)
      integer ii,jj, m, mm, dc(2),df(2,0:rr2-1),iff,jff
      integer refratio

      integer i,j,k, r2,ic,jc, i1, j1


c     # Idea is to fill in a dummy fine grid using zeroth-order
c     # interpolation.  Then assign qc, auxc values

      refratio = 2
      r2 = refratio*refratio

      call build_transform(transform_ptr,a)

c     # This needs to be written for refratios .ne. 2.
      k = 0
      do jj = 0,1
         do ii = 0,1
c           # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

c           # Direction on fine grid (converted using metric). Divide
c           # by refratio to scale length to unit vector
            df(1,k) = (a(1,1)*dc(1) + a(1,2)*dc(2))/refratio
            df(2,k) = (a(2,1)*dc(1) + a(2,2)*dc(2))/refratio

            k = k + 1
         enddo
      enddo
c     # Create map :

      if (idir .eq. 0) then
c        # this ensures that we get 'hanging' corners

         if (iface_coarse .eq. 0) then
            ic = 1
         elseif (iface_coarse .eq. 1) then
            ic = mx
         endif
         do jc = 1,mx
            i1 = ic
            j1 = jc
            call fclaw2d_transform_face_half(i1,j1,i2,j2,
     &            transform_ptr)

            skip_this_grid = .false.
            do k = 0,r2-1
               if (.not. is_valid_correct(idir,i2(k),j2(k),mx,my))
     &               then
                  skip_this_grid = .true.
                  exit
               endif
            enddo
c           # Idea is to fill in ghost cells in groups of four
c           # fine grid ghost cells.
            if (.not. skip_this_grid) then
               do k = 0,rr2-1
                  iff = i2(0) + df(1,k)
                  jff = j2(0) + df(2,k)
                  maskfine(iff,jff) = 1
                  do m = 1,meqn
                     qfine_dummy(iff,jff,m) = qcoarse(ic,jc,m)
                  enddo
                  do m = 1,maux
                     auxfine_dummy(iff,jff,m) = auxcoarse(ic,jc,m)
                  enddo
               enddo
            endif
         enddo
      else
         if (iface_coarse .eq. 2) then
            jc = 1
         elseif (iface_coarse .eq. 3) then
            jc = my
         endif
         do ic = 1,mx
    1       i1 = ic
            j1 = jc
            call fclaw2d_transform_face_half(i1,j1,i2,j2,
     &            transform_ptr)
            skip_this_grid = .false.
            do k = 0,r2-1
               if (.not. is_valid_correct(idir,i2(k),j2(k),mx,my))
     &               then
                  skip_this_grid = .true.
                  exit
               endif
            enddo
            if (.not. skip_this_grid) then
               do k = 0,rr2-1
                  iff = i2(0) + df(1,k)
                  jff = j2(0) + df(2,k)
                  maskfine(iff,jff) = 1
                  do m = 1,meqn
                     qfine_dummy(iff,jff,m) = qcoarse(ic,jc,m)
                  enddo
                  do m = 1,maux
                     auxfine_dummy(iff,jff,m) = auxcoarse(ic,jc,m)
                  enddo
               enddo

            endif                       !! Don't skip this grid
         enddo                          !! i loop
      endif                             !! end idir branch

c     # We only fill in qc arrays which have been set above.
      if (maskfine(0,1) .ne. 0) then
         do j = 1,my
            do mm = 1,meqn
                qc0(j,mm) = qfine_dummy(0,j,mm)
            enddo
            do mm = 1,maux
                auxc0(j,mm) = auxfine_dummy(0,j,mm)
            enddo
         enddo
      elseif (maskfine(mx+1,1) .ne. 0) then
         do j = 1,my
            do mm = 1,meqn
                qc1(j,mm) = qfine_dummy(mx+1,j,mm)
            enddo
            do mm = 1,maux
                auxc1(j,mm) = auxfine_dummy(mx+1,j,mm)
            enddo
         enddo
      elseif (maskfine(1,0) .ne. 0) then
         do i = 1,mx
            do mm = 1,meqn
                qc2(i,mm) = qfine_dummy(i,0,mm)
            enddo
            do mm = 1,maux
                auxc2(i,mm) = auxfine_dummy(i,0,mm)
            enddo
         enddo
      elseif (maskfine(1,my+1) .ne. 0) then
         do i = 1,mx
            do mm = 1,meqn
                qc3(i,mm) = qfine_dummy(i,my+1,mm)
            enddo
            do mm = 1,maux
                auxc3(i,mm) = auxfine_dummy(i,my+1,mm)
            enddo
         enddo
      endif


      end




