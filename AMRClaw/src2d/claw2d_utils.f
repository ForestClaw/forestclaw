      subroutine inputparms(mx_leaf,my_leaf,initial_dt, tfinal,
     &      max_cfl, desired_cfl,nout, src_term, verbose, mcapa,
     &      maux,meqn,mwaves, maxmwaves,mthlim,mbc,mthbc,order,
     &      minlevel,maxlevel,icycle)
      implicit none

      integer mx_leaf, my_leaf
      double precision initial_dt, tfinal, max_cfl, desired_cfl
      integer nout, src_term, mcapa, maux, meqn, mwaves
      integer maxmwaves,mbc, verbose
      integer mthbc(4),mthlim(maxmwaves), order(2)
      integer maxlevel, minlevel, icycle
      logical subcycle

      integer mw, m


      open(55,file='claw2ez.data')

      read(55,*) mx_leaf
      read(55,*) my_leaf

c     timestepping variables
      read(55,*) nout
      read(55,*) tfinal

      read(55,*) initial_dt
      read(55,*) max_cfl
      read(55,*) desired_cfl

      read(55,*) (order(m),m=1,2)
      read(55,*) verbose
      read(55,*) src_term
      read(55,*) mcapa
      read(55,*) maux

      read(55,*) meqn
      read(55,*) mwaves

      if (mwaves .gt. maxmwaves) then
         write(6,*) 'ERROR : (claw_utils.f) mwaves > maxmwaves'
         write(6,*) 'mwaves = ', mwaves
         write(6,*) 'maxmwaves = ', maxmwaves
         stop
      endif

      read(55,*) (mthlim(mw), mw=1,mwaves)

      read(55,*) mbc
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)

      read(55,*) minlevel
      read(55,*) maxlevel

      read(55,*) subcycle
      if (subcycle) then
         icycle = 1
      else
         icycle = 0
      endif


      close(55)
      end


c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine exchange_face_ghost(mx,my,mbc,meqn, qthis,qneighbor,
     &      idir)
      implicit none

      integer mx,my,mbc,meqn,idir
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq

c     # We only have to consider the case of exchanges on
c     # the high side of 'this' grid.
c     # We do need to do a complete exchange, though.
      if (idir .eq. 0) then
         do j = 1,my
            do ibc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # x-direction (idir == 0)
                  qthis(mx+ibc,j,mq) = qneighbor(ibc,j,mq)
                  qneighbor(ibc-mbc,j,mq) = qthis(mx+ibc-mbc,j,mq)
               enddo
            enddo
         enddo
      else
         do i = 1,mx
            do jbc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  qthis(i,my+jbc,mq) = qneighbor(i,jbc,mq)
                  qneighbor(i,jbc-mbc,mq) = qthis(i,my-mbc+jbc,mq)
               enddo
            enddo
         enddo
      endif
      end



c     # Average fine grid to a coarse grid neighbor or copy from neighboring
c     # grid at same level.
      subroutine average_face_ghost(mx,my,mbc,meqn,
     &      qfine,qcoarse,idir,iside,num_neighbors,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iside, num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer mq,r2
      integer i, i1, i2, ibc, ii, ii1
      integer j, j1, j2, jbc, jj, jj1

      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
      if (idir .eq. 0) then
         j1 = igrid*my/num_neighbors+1
         j2 = (igrid+1)*my/num_neighbors
         do j = j1, j2
            do ibc = 1,mbc
               do mq = 1,meqn
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ii1 = (ibc-1)*refratio + ii
                        jj1 = (j-1)*refratio + jj
                        if (iside .eq. 0) then
                           sum = sum + qfine(mx-ii1+1,jj1,mq)
                        else
                           sum = sum + qfine(ii1,jj1,mq)
                        endif
                     enddo
                  enddo
                  if (iside .eq. 0) then
                     qcoarse(1-ibc,j,mq) = sum/r2
                  else
                     qcoarse(mx+ibc,j,mq) = sum/r2
                  endif
               enddo
            enddo
         enddo
      else
         i1 = igrid*mx/num_neighbors+1
         i2 = (igrid+1)*mx/refratio
         do i = i1,i2
            do jbc = 1,mbc
               do mq = 1,meqn
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ii1 = (i-1)*refratio + ii
                        jj1 = (jbc-1)*refratio + jj
                        if (iside .eq. 2) then
                           sum = sum + qfine(ii1,my-jj1+1,mq)
                        else
                           sum = sum + qfine(ii1,jj1,mq)
                        endif
                     enddo
                  enddo
                  if (iside .eq. 0) then
                     qcoarse(i,1-jbc,mq) = sum/r2
                  else
                     qcoarse(i,my+jbc,mq) = sum/r2
                  endif
               enddo
            enddo
         enddo
      endif

      end

      subroutine interpolate_face_ghost(mx,my,mbc,meqn,qfine,qcoarse,
     &      idir,iside,num_neighbors,refratio,igrid)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iside,num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2
      integer i, i1, i2, ibc, ii, ii1
      integer j, j1, j2, jbc, jj, jj1


c     # Assign values to fine grid ghost.  For now, just copy coarse grid
c     values; Get linear part later.
      if (idir .eq. 0) then
         j1 = igrid*my/num_neighbors+1
         j2 = (igrid+1)*my/num_neighbors
         do j = j1, j2
            do ibc = 1,mbc
               do mq = 1,meqn
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ii1 = (ibc-1)*refratio + ii
                        jj1 = (j-1)*refratio + jj
                        if (iside .eq. 0) then
                           qfine(ii1-mbc,jj1,mq) = qcoarse(ibc,j,mq)
                        elseif (iside .eq. 1) then
                           qfine(mx+ii1,jj1,mq) = qcoarse(mx+ibc,j,mq)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         i1 = igrid*mx/num_neighbors+1
         i2 = (igrid+1)*mx/refratio
         do i = i1,i2
            do jbc = 1,mbc
               do mq = 1,meqn
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ii1 = (i-1)*refratio + ii
                        jj1 = (jbc-1)*refratio + jj
                        if (iside .eq. 2) then
                           qfine(ii1,jj1-mbc,mq) = qcoarse(i,jbc,mq)
                        elseif (iside .eq. 3) then
                           qfine(ii1,my+jj1,mq) = qcoarse(i,my+jbc,mq)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif


      end


c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_corner(mx,my,mbc,meqn,
     &      qcoarse,qfine,icorner,refratio)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2

      r2 = refratio*refratio

      if (refratio .eq. 1) then
c        # We only need to worry about corners 1 and 3 (lr and ur).
c        # for the  complete exchange.  The other corners are somebody
c        # else's (lr,ur) corners.
         do mq = 1,meqn
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Exchange corner information at boundaries.
                  if (icorner .eq. 1) then
c                    # Fix this!
                     qcoarse(mx+ibc,jbc-mbc,mq) =
     &                     qfine(ibc,my+jbc-mbc,mq)
                     qfine(ibc-mbc,my+jbc,mq) =
     &                     qcoarse(mx+ibc-mbc,jbc,mq)
                  elseif (icorner .eq. 3) then
                     qcoarse(mx+ibc,my+jbc,mq) =
     &                     qfine(ibc,jbc,mq)
                     qfine(ibc-mbc,jbc-mbc,mq) =
     &                     qcoarse(mx+ibc-mbc,my+jbc-mbc,mq)
                  endif
               enddo
            enddo
         enddo
      else
c        # Average fine grid onto coarse grid
         write(6,'(A,A)') 'average_corner_step1 : fine grid ',
     &         ' averaging at corners not yet implemented'
         stop
         do mq = 1,meqn
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (icorner .eq. 0) then
                  elseif (icorner .eq. 1) then
                  elseif (icorner .eq. 2) then
                  else
                  endif
               enddo
            enddo
         enddo
      endif

      end


      subroutine set_phys_corner_ghost(mx,my,mbc,meqn,q,icorner,t,dt,
     &      mthbc)
      implicit none

      integer mx,my,mbc,meqn,icorner, mthbc(4)
      double precision t,dt
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

c     # Do something here....

      end


      subroutine exchange_phys_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iside)
      implicit none

      integer mx, my, mbc, meqn, iside, icorner
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ibc, jbc, mq

c     # Fill in corner ghost cells that overlap the physical boundary. In this
c     case, the corner ghost cells are copied from a face neighbor.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (iside .eq. 1) then
                  if (icorner .eq. 1) then
                     qthis(mx+ibc,jbc-mbc,mq) =
     &                     qneighbor(ibc,jbc-mbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(mx+ibc-mbc,jbc-mbc,mq)
                  elseif(icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(ibc,my+jbc,mq)
                     qneighbor(ibc-mbc,my+jbc,mq) =
     &                     qthis(mx+ibc-mbc,my+jbc,mq)
                  endif
               elseif (iside .eq. 3) then
                  if (icorner .eq. 2) then
                     qthis(ibc-mbc,my+jbc,mq) =
     &                     qneighbor(ibc-mbc,jbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(ibc-mbc,my+jbc-mbc,mq)
                  elseif(icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(mx+ibc,jbc,mq)
                     qneighbor(mx+ibc,jbc-mbc,mq) =
     &                     qthis(mx+ibc,my+jbc-mbc,mq)
                  endif
               endif
            enddo
         enddo
      enddo



      end


      subroutine exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner)
      implicit none

      integer mx, my, mbc, meqn, icorner
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

c     # NOTE : qneighbor is not yet a valid pointer.

c     # Only exchanging high side corners

c     # We only need to worry about corners 1 and 3 (lr and ur).
c     # for the  complete exchange.  The other corners are somebody
c     # else's (lr,ur) corners.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Exchange corner information at boundaries.
c              # Do this until we get the corner ids fixed
               if (icorner .eq. 1) then
                  qthis(mx+ibc,jbc-mbc,mq) =
     &                  qneighbor(ibc,my+jbc-mbc,mq)
                  qneighbor(ibc-mbc,my+jbc,mq) =
     &                  qthis(mx+ibc-mbc,jbc,mq)
               elseif (icorner .eq. 3) then
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(ibc,jbc,mq)
                  qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                  qthis(mx+ibc-mbc,my+jbc-mbc,mq)
               endif
            enddo
         enddo
      enddo


      end


      subroutine compute_sum(mx,my,mbc,meqn,dx,dy,q,sum)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, sum
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j

      sum = 0
      do i = 1,mx
         do j = 1,my
            sum = sum + q(i,j,1)*dx*dy
         enddo
      enddo
      end


c     # Conservative intepolation to fine grid patch
      subroutine interpolate_to_fine_patch(mx,my,mbc,meqn, qcoarse,
     &      qfine, p4est_refineFactor,refratio,fine_patch_idx)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer fine_patch_idx
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      double precision qc, shiftx, shifty, sl, sr, gradx, grady
      double precision compute_slopes


c     # Use limiting done in AMRClaw.
      mth = 5

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(fine_patch_idx,refratio)
      jg = (fine_patch_idx-ig)/refratio

      i1 = ig*mx/p4est_refineFactor + 1
      i2 = (ig+1)*mx/p4est_refineFactor

      j1 = jg*my/p4est_refineFactor + 1
      j2 = (jg+1)*my/p4est_refineFactor

      do mq = 1,meqn
         do i = i1,i2
            do j = j1,j2
               qc = qcoarse(i,j,mq)

c              # Compute limited slopes in both x and y. Note we are not
c              # really computing slopes, but rather just differences.
c              # Scaling is accounted for in 'shiftx' and 'shifty', below.
               sl = (qc - qcoarse(i-1,j,mq))
               sr = (qcoarse(i+1,j,mq) - qc)
               gradx = compute_slopes(sl,sr,mth)
               sl = (qc - qcoarse(i,j-1,mq))
               sr = (qcoarse(i,j+1,mq) - qc)
               grady = compute_slopes(sl,sr,mth)

c              # Fill in refined values on coarse grid cell (i,j)
               do ii = 1,refratio
                  do jj = 1,refratio
                     shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                     shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                     qfine((i-1)*refratio + ii,(j-1)*refratio + jj,mq)
     &                     = qc + shiftx*gradx + shifty*grady
                  enddo
               enddo
            enddo
         enddo
      enddo

      end

      double precision function compute_slopes(sl,sr,mth)
      implicit none

      double precision sl,sr, s, sc, limiter, slim
      integer mth

c     # ------------------------------------------------
c     # Slope limiting done in amrclaw - see filpatch.f
c       dupc = valp10 - valc
c       dumc = valc   - valm10
c       ducc = valp10 - valm10
c       du   = dmin1(dabs(dupc),dabs(dumc))        <--
c       du   = dmin1(2.d0*du,.5d0*dabs(ducc))      <-- Not quite sure I follow
c
c       fu = dmax1(0.d0,dsign(1.d0,dupc*dumc))
c
c       dvpc = valp01 - valc
c       dvmc = valc   - valm01
c       dvcc = valp01 - valm01
c       dv   = dmin1(dabs(dvpc),dabs(dvmc))
c       dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
c       fv = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))
c
c       valint = valc + eta1*du*dsign(1.d0,ducc)*fu
c      .      + eta2*dv*dsign(1.d0,dvcc)*fv
c     # ------------------------------------------------

c     # ------------------------------------------------
c     # To see what Chombo does, look in InterpF.ChF
c     # (in Chombo/lib/src/AMRTools) for routine 'interplimit'
c     # Good luck.
c     # ------------------------------------------------

      if (mth .le. 4) then
c        # Use minmod, superbee, etc.
         slim = limiter(sl,sr,mth)
         compute_slopes = slim*sl
      else
c        # Use AMRClaw slopes
         sc = (sl + sr)/2.d0
         compute_slopes = min(abs(sl),abs(sr),abs(sc))*
     &         max(0.d0,sign(1.d0,sl*sr))*sign(1.d0,sc)
      endif


      end
