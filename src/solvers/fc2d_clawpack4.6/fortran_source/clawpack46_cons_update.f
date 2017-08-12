c    # -----------------------------------------------------------------     
c    # These are called from ghost filling routines
c    # -----------------------------------------------------------------     


c    # Called from averaging routine, using fine grid (possibly from 
c    # a ghost patch) to compute correction terms.  
      subroutine clawpack46_fort_cons_coarse_correct(mx,my,mbc,meqn,
     &                                          idir, iface_coarse,
     &                                          area0, area1, 
     &                                          area2, area3,
     &                                          qcoarse,
     &                                          fpcoarse0,fmcoarse1,
     &                                          gpcoarse2,gmcoarse3,
     &                                          fm0, fp1,gm0, gp1,
     &                                          rp0, rp1, rp2, rp3,
     &                                          maskfine, qfine_dummy,
     &                                          transform_cptr)

      implicit none

      integer mx,my,mbc,meqn,idir,iface_coarse

      integer maskfine(1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine_dummy(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision area0(my), area1(my), area2(mx), area3(mx)

      double precision fpcoarse0(my,meqn)
      double precision fmcoarse1(my,meqn)
      double precision gpcoarse2(mx,meqn)
      double precision gmcoarse3(mx,meqn)

      double precision fm0(my,meqn)
      double precision fp1(my,meqn)
      double precision gm0(mx,meqn)
      double precision gp1(mx,meqn)

      double precision rp0(my,meqn)
      double precision rp1(my,meqn)
      double precision rp2(mx,meqn)
      double precision rp3(mx,meqn)

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


c     # Visit each fine grid face and fill in correction term into
c     # qfine_dummy. We don't really know which face corresponds
c     # to the common face shared between fine and coarse grids, so we
c     # do them all.

      do mq = 1,meqn
         do j = 1,my/2
            jj1 = 2*(j-1) + 1     !! indexing starts at 1
            jj2 = 2*(j-1) + 2
            deltam = (fm0(jj1,mq) + fm0(jj2,mq)) +
     &               (rp0(jj1,mq) + rp0(jj2,mq))
            deltap = (fp1(jj1,mq) + fp1(jj2,mq)) +
     &               (rp1(jj1,mq) + rp1(jj2,mq))

            do jj = jj1,jj2
               do ii = 1,mbc
                  qfine_dummy(1-ii,jj,mq) = deltam
                  qfine_dummy(mx+ii,jj,mq) = deltap
               enddo
            enddo
         enddo

         do i = 1,mx/2
            ii1 = 2*(i-1) + 1        !! indexing starts at 1
            ii2 = 2*(i-1) + 2
            deltam = (gm0(ii1,mq) + gm0(ii2,mq) +
     &               rp2(ii1,mq) + rp2(ii2,mq))
            deltap = (gp1(ii1,mq) + gp1(ii2,mq) +
     &               rp3(ii1,mq) + rp3(ii2,mq))
            do ii = ii1,ii2
               do jj = 1,mbc
                  qfine_dummy(ii,1-jj,mq) = deltam
                  qfine_dummy(ii,my+jj,mq) = deltap
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
                     deltac = fineval - fpcoarse0(jc,mq) 
                     areac = area0(jc)
                  else
                     deltac = fineval - fmcoarse1(jc,mq)
                     areac = area1(jc)
                  endif
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) - deltac/areac

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
                     deltac = fineval - gpcoarse2(ic,mq)
                     areac = area2(ic)
                  else               
                     deltac = fineval - gmcoarse3(ic,mq) 
                     areac = area3(ic)
                  endif
                  qcoarse(ic,jc,mq) = qcoarse(ic,jc,mq) - deltac/areac
               endif                    !! skip grid loop
            enddo
         enddo
      endif


      end
   
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


c    # -----------------------------------------------------------------     
c    # These are called from step2 routine (after update)
c    # -----------------------------------------------------------------     

      subroutine clawpack46_accumulate_cons_updates(mx,my,mbc,meqn,
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


      subroutine clawpack46_accumulate_riemann_problem(mx,my,mbc,meqn,
     &      maux, dt,
     &      el0, el1, el2, el3,
     &      mside,idir,iside,qc,auxc,qfine,auxfine,rp_accum,
     &      rpn2_cons,ql,qr,auxl,auxr,flux_diff)

      implicit none

      integer mx,my,mbc,meqn, maux, iside, mside, idir
      double precision dt

      double precision el0(my), el1(my), el2(mx), el3(mx)

      double precision qc(mside,meqn)
      double precision auxc(mside,maux)
      double precision rp_accum(mside,meqn)

      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision ql(meqn),qr(meqn),auxl(maux), auxr(maux)
      double precision flux_diff(meqn)

      integer i,j,m
      double precision elen

      if (idir .eq. 0) then
         do j = 1,my
            if (iside .eq. 0) then
c              # left side
               do m = 1,maux
                  auxl(m) = auxc(j,m)
                  auxr(m) = auxfine(0,j,m)
               enddo
               do m = 1,meqn
                  ql(m) = qc(j,m)
                  qr(m) = qfine(0,j,m)
               enddo
               elen = el0(j)
            elseif (iside .eq. 1) then
c              # right side
               do m = 1,maux
                  auxl(m) = auxfine(mx+1,j,m)
                  auxr(m) = auxc(j,m)
               enddo
               do m = 1,meqn
                  ql(m) = qfine(mx+1,j,m)
                  qr(m) = qc(j,m)
               enddo
               elen = el1(j)
            endif

            call rpn2_cons(meqn,maux,idir,ql,qr,auxl,auxr,flux_diff)

            do m = 1,meqn
               rp_accum(j,m) = rp_accum(j,m) + dt*elen*flux_diff(m)
            enddo
         enddo
      elseif (idir .eq. 1) then
         do i = 1,mx
            if (iside .eq. 2) then
c              # bottom side
               do m = 1,maux
                  auxl(m) = auxc(i,m)
                  auxr(m) = auxfine(i,0,m)
               enddo
               do m = 1,meqn
                  ql(m) = qc(i,m)
                  qr(m) = qfine(i,0,m)
               enddo
               elen = el2(i)
            elseif (iside .eq. 3) then
c              # top side
               do m = 1,maux
                  auxl(m) = auxfine(i,my+1,m)
                  auxr(m) = auxc(i,m)
               enddo
               do m = 1,meqn
                  ql(m) = qfine(i,my+1,m)
                  qr(m) = qc(i,m)
               enddo
               elen = el3(i)
            endif

            call rpn2_cons(meqn,maux,idir,ql,qr,auxl,auxr,flux_diff)

            do m = 1,meqn
               rp_accum(i,m) = rp_accum(i,m) + dt*elen*flux_diff(m)
            enddo
         enddo
      endif

      end



      logical function is_valid_correct(idir,i,j,mx,my)
      implicit none

      integer i,j,mx,my, idir
      logical i1, j1

      i1 = 1 .le. i .and. i .le. mx
      j1 = 1 .le. j .and. j .le. my

      is_valid_correct = xor(i1,j1)


      end
