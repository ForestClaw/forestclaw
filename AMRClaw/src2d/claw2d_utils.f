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

c    # ----------------------------------------------------------------------------------
c    # Internal boundary conditions
c    # ----------------------------------------------------------------------------------

c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine exchange_face_ghost(mx,my,mbc,meqn,qthis,qneighbor,
     &      idir)
      implicit none

      integer mx,my,mbc,meqn,idir
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq

c     # High side of 'qthis' exchanges with low side of
c     # 'qneighbor'
      if (idir .eq. 0) then
         do j = 1,my
            do ibc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # x-direction (idir == 0)
                  qthis(mx+ibc,j,mq) = qneighbor(ibc,j,mq)
                  qneighbor(1-ibc,j,mq) = qthis(mx+1-ibc,j,mq)
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


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iside'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost(mx,my,mbc,meqn,qcoarse,
     &      qfine,idir,iface_coarse,num_neighbors,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer mq,r2
      integer i, i1, ibc, ii, ifine
      integer j, j1, jbc, jj, jfine

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
      if (idir .eq. 0) then
         j1 = igrid*my/num_neighbors
         do mq = 1,meqn
            do j = 1,my/num_neighbors
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (j-1)*refratio + jj
                        if (iface_coarse .eq. 0) then
                           sum = sum + qfine(mx-ifine+1,jfine,mq)
                        else
                           sum = sum + qfine(ifine,jfine,mq)
                        endif
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     qcoarse(1-ibc,j+j1,mq) = sum/r2
                  else
                     qcoarse(mx+ibc,j+j1,mq) = sum/r2
                  endif
               enddo
            enddo
         enddo
      else
         i1 = igrid*mx/num_neighbors
         do i = 1,mx/num_neighbors
            do jbc = 1,mbc
               do mq = 1,meqn
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (i-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (iface_coarse .eq. 2) then
                           sum = sum + qfine(ifine,my-jfine+1,mq)
                        else
                           sum = sum + qfine(ifine,jfine,mq)
                        endif
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     qcoarse(i+i1,1-jbc,mq) = sum/r2
                  else
                     qcoarse(i+i1,my+jbc,mq) = sum/r2
                  endif
               enddo
            enddo
         enddo
      endif

      end

      subroutine interpolate_face_ghost(mx,my,mbc,meqn,qcoarse,qfine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2
      integer i, i1, i2, ibc, ii, ii1, ifine
      integer j, j1, j2, jbc, jj, jj1, jfine
      integer ic_add, jc_add
      integer linear_terms

c     # To be figured out later
      linear_terms = 0

c     # 'iface_coarse is relative to the coarse grid


c     # Assign values to fine grid ghost.  For now, just copy coarse grid
c     # values; Get linear part later.
      do mq = 1,meqn
         if (idir .eq. 0) then
            jc_add = igrid*my/num_neighbors
c           # Need to also get corners that are at the center of a coarse
c           # face
            do j = (1-igrid),my/num_neighbors + (1-igrid)*mbc
               do ibc = 1,mbc
                  do jj = 1,refratio
                     ifine = ibc
                     jfine = (j-1)*refratio + jj
                     if (iface_coarse .eq. 0) then
c                       # qfine is at left edge of coarse grid
                        qfine(mx+ifine,jfine,mq) =
     &                        qcoarse(1,j+jc_add,mq) + linear_terms
                     elseif (iface_coarse .eq. 1) then
c                       # qfine is at right edge of coarse grid
                        qfine(1-ifine,jfine,mq) =
     &                        qcoarse(mx,j+jc_add,mq) + linear_terms
                     endif
                  enddo
               enddo
            enddo
         else
            ic_add = igrid*mx/num_neighbors
            do i =(1-igrid),mx/num_neighbors + (1-igrid)*mbc
               do ii = 1,refratio
                  do jbc = 1,mbc
                     ifine = (i-1)*refratio + ii
                     jfine = jbc
                     if (iface_coarse .eq. 2) then
c                       # qfine is at bottom edge of coarse grid
                        qfine(ifine,my+jfine,mq) =
     &                        qcoarse(i+ic_add,1,mq) + linear_terms
                     else if (iface_coarse .eq. 3) then
c                       # qfine is at top edge of coarse grid
                        qfine(ifine,1-jfine,mq) =
     &                        qcoarse(i+ic_add,my,mq) + linear_terms
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
      end

      subroutine exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner_this)
      implicit none

      integer mx, my, mbc, meqn, icorner_this
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

c     # Only exchanging high side corners

c     # We only need to worry about corners 1 and 3 (lr and ur).
c     # for the  complete exchange.  The other corners are somebody
c     # else's (lr,ur) corners.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Exchange initiated only at high side (1,3) corners
               if (icorner_this .eq. 1) then
                  qthis(mx+ibc,1-jbc,mq) =
     &                  qneighbor(ibc,my+1-jbc,mq)
                  qneighbor(1-ibc,my+jbc,mq) =
     &                  qthis(mx+1-ibc,jbc,mq)
               elseif (icorner_this .eq. 3) then
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(ibc,jbc,mq)
                  qneighbor(1-ibc,1-jbc,mq) =
     &                  qthis(mx+1-ibc,my+1-jbc,mq)
               endif
            enddo
         enddo
      enddo
      end



c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_corner_ghost(mx,my,mbc,meqn,
     &      refratio,qcoarse,qfine,icorner_coarse)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner_coarse
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2
      integer ifine, jfine

      r2 = refratio*refratio
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Average fine grid corners onto coarse grid ghost corners
               sum = 0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (ibc-1)*refratio + ii
                     jfine = (jbc-1)*refratio + jj
                     if (icorner_coarse .eq. 0) then
                        sum = sum + qfine(mx+1-ifine,my+1-jfine,mq)
                     elseif (icorner_coarse .eq. 1) then
                        sum = sum + qfine(ifine,my+1-jfine,mq)
                     elseif (icorner_coarse .eq. 2) then
                        sum = sum + qfine(mx+1-ifine,jfine,mq)
                     elseif (icorner_coarse .eq. 3) then
                        sum = sum + qfine(ifine,jfine,mq)
                     endif
                  enddo
               enddo
               if (icorner_coarse .eq. 0) then
                  qcoarse(1-ibc,1-jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 1) then
                  qcoarse(mx+ibc,1-jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 2) then
                  qcoarse(1-ibc,my+jbc,mq) = sum/r2
               elseif (icorner_coarse .eq. 3) then
                  qcoarse(mx+ibc,my+jbc,mq) = sum/r2
               endif
            enddo
         enddo
      enddo

      end

      subroutine interpolate_corner_ghost(mx,my,mbc,meqn,refratio,
     &      qcoarse,qfine,icorner_coarse)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ic, jc, mq, ibc,jbc, mth,i,j
      double precision qc, sl, sr, gradx, grady, shiftx, shifty
      double precision compute_slopes, value

      mth = 5

      if (icorner_coarse .eq. 0) then
         ic = 1
         jc = 1
      elseif (icorner_coarse .eq. 1) then
         ic = mx
         jc = 1
      elseif (icorner_coarse .eq. 2) then
         ic = 1
         jc = my
      elseif (icorner_coarse .eq. 3) then
         ic = mx
         jc = my
      endif

      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)
c        # Interpolate coarse grid corners to fine grid corner ghost cells
         do ibc = 1,mbc
            do jbc = 1,mbc

c              # Compute limited slopes in both x and y. Note we are not
c              # really computing slopes, but rather just differences.
c              # Scaling is accounted for in 'shiftx' and 'shifty', below.
               sl = (qc - qcoarse(ic-1,jc,mq))
               sr = (qcoarse(ic+1,jc,mq) - qc)
               gradx = compute_slopes(sl,sr,mth)

               sl = (qc - qcoarse(ic,jc-1,mq))
               sr = (qcoarse(ic,jc+1,mq) - qc)
               grady = compute_slopes(sl,sr,mth)

c              # Fill in interpolated values on fine grid cell
               shiftx = (ibc - refratio/2.d0 - 0.5d0)/refratio
               shifty = (jbc - refratio/2.d0 - 0.5d0)/refratio

               value = (qc + shiftx*gradx + shifty*grady)
               if (icorner_coarse .eq. 0) then
                  qfine(mx+ibc,my+jbc,mq) = value
               elseif (icorner_coarse .eq. 1) then
                  qfine(mx+ibc,1-jbc,mq) = value
               elseif (icorner_coarse .eq. 2) then
                  qfine(1-ibc,my+jbc,mq) = value
               else
                  qfine(1-ibc,1-jbc,mq) = value
               endif
            enddo
         enddo
      enddo

      end

c     # ----------------------------------------------------------------------
c     # Physical boundary conditions
c     # ----------------------------------------------------------------------
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
                  elseif (icorner .eq. 3) then
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


c     # Use bc2 to assign boundary conditions to corners with the physical
c     # boundary of another patch.   Use this routine instead of the exchange,
c     above, in the case where q and its neighbor are not at the same level.
c     In this case, we can't just copy, and intepolation/averaging doesn't
c     make sense, since we have a physical boundary condition
      subroutine set_phys_interior_corner_ghost(mx,my,mbc,meqn, q,
     &      icorner, iface, mthbc)
      implicit none

      integer mx, my, mbc, meqn, icorner, iface, mthbc
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer this_side_bc, mthbc_tmp

c     # Redo this side, knowing that it now has valid interior ghost cell data
c     # from which we can copy.

      end

c    # ---------------------------------------------------------------------
c    # Tagging for refinement/coarsening
c    # ---------------------------------------------------------------------

      subroutine tag_for_refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      q,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      double precision xlower, ylower, dx, dy
      double precision q(1-mx:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j
      double precision xc,yc

      tag_patch = 0
      do i = 1,mx
         do j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            if (abs(xc - 0.5d0) < dx) then
               tag_patch = 1
            endif
         enddo
      enddo

      end

c     # Conservative intepolation to fine grid patch
      subroutine interpolate_to_fine_patch(mx,my,mbc,meqn, qcoarse,
     &      qfine, p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
      double precision qc, shiftx, shifty, sl, sr, gradx, grady
      double precision compute_slopes


c     # Use limiting done in AMRClaw.
      mth = 5

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

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

      subroutine average_to_coarse_patch(mx,my,mbc,meqn,qcoarse,qfine,
     &      p4est_refineFactor, refratio, igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, ig, jg, i1, i2, j1, j2, ii, jj, ifine, jfine
      integer mq
      double precision sum, r2

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      i1 = ig*mx/p4est_refineFactor + 1
      i2 = (ig+1)*mx/p4est_refineFactor

      j1 = jg*my/p4est_refineFactor + 1
      j2 = (jg+1)*my/p4est_refineFactor

      r2 = refratio*refratio
      do mq = 1,meqn
         do j = j1,j2
            do i = i1,i2
               sum = 0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     sum = sum + qfine(ifine,jfine,mq)
                  enddo
               enddo
               qcoarse(i,j,mq) = sum/r2
            enddo
         enddo
      enddo
      end


c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
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
