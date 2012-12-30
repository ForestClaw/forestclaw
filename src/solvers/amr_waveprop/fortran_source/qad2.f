      subroutine qad2(maxm,meqn,mwaves,mbc,mx,my,maux,
     &      qold,qold_coarse,auxold, auxold_coarse,intersectsBoundary,
     &      mxc,myc,dt,dx,dy, ql,qr,auxl,auxr, wave,s,
     &      amdq,apdq, qadd_x,qadd_y,auxtype)


      implicit none
c     # Input parameters
      integer maxm,meqn,mwaves,mbc,mx,my,mxc,myc,maux
      integer intersectsBoundary(4)
      double precision dt, dx, dy

c     # note that mbc is only relevant for the fine grid;  the coarse
c     # grid only has a single layer of ghost cells.
      double precision qold(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qold_coarse(0:mxc+1,0:myc+1,meqn)
      double precision auxold(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision auxold_coarse(0:mxc+1,0:myc+1,maux)

      double precision qadd_x(mx+1,my,meqn)
      double precision qadd_y(mx,my+1,meqn)

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)

      double precision auxl(1-mbc:maxm+mbc,maux)
      double precision auxr(1-mbc:maxm+mbc,maux)
      logical inqad

      common /comqad/ inqad

      character*10 auxtype(maux)

      integer i,j,ic,jc,ifine,jfine,ir,jr,m, iref_ratio, ixy
      logical ib(4)

      integer iunit

      inqad = .true.


      iref_ratio = mx/mxc
      if (iref_ratio == 1) then
c        # Should be equvilient to
c        #           (ib(1) .and. ib(2) .and. ib(3) .and. ib(4))
c        #
c        # We are on the coarsest grid,and there is really nothing we
c        # should be doing here, except perhaps zeroing out qadd_coarse
c        # which we did above.
         inqad = .false.
         return
      endif


c     # Now check boundaries to get something reasonable if the
c     # coarse grid intersects a physical boundary.
c     # Note that it doesn't really matter what we put here, but
c     # we want something that makes sense for the problem so that
c     # when we pass everything to the Riemann solver, we don't
c     # crash with non-sensical values.

      do i = 1,4
         ib(i) = intersectsBoundary(i) == 1
      enddo

      if (ib(1)) then
         do  jc = 0,myc+1
            do m = 1,meqn
               qold_coarse(0,jc,m) = qold_coarse(1,jc,m)
            enddo
            do m = 1,maux
               auxold_coarse(0,jc,m) = auxold_coarse(1,jc,m)
            enddo
         enddo
      endif

      if (ib(2)) then
         do  jc = 0,myc+1
            do m = 1,meqn
               qold_coarse(mxc+1,jc,m) = qold_coarse(mxc,jc,m)
            enddo
            do m = 1,maux
               auxold_coarse(mxc+1,jc,m) = auxold_coarse(mxc,jc,m)
            enddo
         enddo
      endif

      if (ib(3)) then
         do  ic = 0,mxc+1
            do m = 1,meqn
               qold_coarse(ic,0,m) = qold_coarse(ic,1,m)
            enddo
            do m = 1,maux
               auxold_coarse(ic,0,m) = auxold_coarse(ic,1,m)
            enddo
         enddo
      endif

      if (ib(4)) then
         do  ic = 0,mxc+1
            do m = 1,meqn
               qold_coarse(ic,myc+1,m) = qold_coarse(ic,myc,m)
            enddo
            do m = 1,maux
               auxold_coarse(ic,myc+1,m) = auxold_coarse(ic,myc,m)
            enddo
         enddo
      endif


c ----------------------------------------------------------------
c     # Side 1
c     # Left edge (xleft)
c ----------------------------------------------------------------

      if (.not. ib(1)) then
         ic = 0
         ifine = 0
c        # Do not loop 1-mbc, myc+mbc - the coarse grid only has
c        # a single layer of ghost cells.
c
c        # We loop over ghost cell values here because the Riemann
c        # solver will loop over these cells and we don't want it to
c        # crash with bogus values.
         do jc = 0,myc+1
c           # Loop over number of cells in fine grid per coarse grid cell
            do jr = 1,iref_ratio
c              # maux might be 0, in which case this loop is not executed.
               jfine = (jc - 1)*iref_ratio + jr
               if (jfine .lt. -1 .or. jfine .gt. my+mbc) cycle
               do m = 1,maux
                  if (auxtype(m) .eq. "xleft") then
                     auxl(jfine,m) = auxold(ifine+1,jfine,m)
                  else
                     auxl(jfine,m) = auxold(ifine,jfine,m)
                  endif
                  if (jfine .gt. -1) then
                     auxr(jfine-1,m) = auxold_coarse(ic,jc,m)
                  endif
               enddo
               do m = 1,meqn
                  ql(jfine,m) = qold(ifine,jfine,m)
                  if (jfine .gt. -1) then
                     qr(jfine-1,m) = qold_coarse(ic,jc,m)
                  endif
               enddo
            enddo
         enddo

c        # Doing an x sweep
         ixy = 1
         call rpn2(ixy,maxm,meqn,mwaves,mbc,my,ql,qr,auxl,auxr,
     &         wave,s,amdq,apdq)

         do m = 1,meqn
            do jc = 1,myc
               do jr = 1,iref_ratio
                  jfine = (jc-1)*iref_ratio + jr
                  if (jfine .lt. -1 .or. jfine .gt. my+mbc) cycle
                  qadd_x(1,jfine,m) = (amdq(jfine,m) + apdq(jfine,m))
               enddo
            enddo
         enddo
      endif
c ----------------------------------------------------------------
c     # Side 2
c     # left edge (xright)
c ----------------------------------------------------------------
      if (.not. ib(2)) then
         ic = mxc + 1
         ifine = mx + 1
c        # Do not loop 1-mbc, myc+mbc - the coarse grid only has
c        # a single layter of ghost cells.
c        # See note for side 1 as to why we loop over ghost cells here...
         do jc = 0,myc+1
c           # Avoid sending bogus values to Riemann solver, if possible.
c           # Loop over number of cells in fine grid per coarse grid cell
            do jr = 1,iref_ratio
               jfine = (jc - 1)*iref_ratio + jr
               if (jfine .lt. -1 .or. jfine .gt. my+mbc) cycle
               do m = 1,maux
                  auxl(jfine,m) = auxold_coarse(ic,jc,m)
                  if (jfine .gt. -1) then
                     auxr(jfine-1,m) = auxold(ifine,jfine,m)
                  endif
               enddo
               do m = 1,meqn
                  ql(jfine,m) = qold_coarse(ic,jc,m)
                  if (jfine .gt. -1) then
                     qr(jfine-1,m) = qold(ifine,jfine,m)
                  endif
               enddo
            enddo
         enddo

c        # Doing an x-sweep
         ixy = 1
         call rpn2(ixy,maxm,meqn,mwaves,mbc,my,ql,qr,auxl,auxr,
     &         wave,s,amdq,apdq)

         do m = 1,meqn
            do jc = 1,myc
               do jr = 1,iref_ratio
                  jfine = (jc-1)*iref_ratio + jr
                  if (jfine .lt. -1 .or. jfine .gt. my+mbc) cycle
                  qadd_x(mx+1,jfine,m) = (amdq(jfine,m) + apdq(jfine,m))
               enddo
            enddo
         enddo
      endif
c ----------------------------------------------------------------
c     # Side 3
c     # Lower edge (yleft)
c ----------------------------------------------------------------
      if (.not. ib(3)) then
         jc = 0
         jfine = 0
c        # Do not loop 1-mbc, mxc+mbc - the coarse grid only has
c        # a single layter of ghost cells.
c        # See note for side 1 as to why we loop over ghost cells here...
         do ic = 0,mxc+1
c           # Avoid sending bogus values to Riemann solver, if possible.
c           # These values aren't really right, but they won't get used.
            do ir = 1,iref_ratio
               ifine = (ic-1)*iref_ratio + ir
               if (ifine .lt. -1 .or. ifine .gt. mx+mbc) cycle
               do m = 1,maux
                  if (auxtype(m) .eq. "yleft") then
                     auxl(ifine,m) = auxold(ifine,jfine+1,m)
                  else
                     auxl(ifine,m) = auxold(ifine,jfine,m)
                  endif
                  if (ifine .gt. -1) then
                     auxr(ifine-1,m) = auxold_coarse(ic,jc,m)
                  endif
               enddo
               do m = 1,meqn
                  ql(ifine,m) = qold(ifine,jfine,m)
                  if (ifine .gt. -1) then
                     qr(ifine-1,m) = qold_coarse(ic,jc,m)
                  endif
               enddo
            enddo
         enddo

         ixy = 2
         call rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &         wave,s,amdq,apdq)

         do m = 1,meqn
            do ic = 1,mxc
               do ir = 1,iref_ratio
                  ifine = (ic-1)*iref_ratio + ir
                  if (ifine .lt. -1 .or. ifine .gt. mx+mbc) cycle
                  qadd_y(ifine,1,m) = (amdq(ifine,m) + apdq(ifine,m))
               enddo
            enddo
         enddo
      endif

c ----------------------------------------------------------------
c     # Side 4
c     # Upper edge (yright)
c ----------------------------------------------------------------
      if (.not. ib(4)) then
         jc = myc+1
         jfine = my+1
c        # Do not loop 1-mbc, mxc+mbc - the coarse grid only has
c        # a single layter of ghost cells.
c        # See note for side 1 as to why we loop over ghost cells here...
         do ic = 0,mxc+1
c           # Avoid sending bogus values to Riemann solver, if possible.
c           # These values aren't really correct, but they won't get used.
            do ir = 1,iref_ratio
               ifine = (ic-1)*iref_ratio + ir
               if (ifine .lt. -1 .or. ifine .gt. mx+mbc) cycle
               do m = 1,maux
                  auxl(ifine,m) = auxold_coarse(ic,jc,m)
                  if (ifine .gt. -1) then
                     auxr(ifine-1,m) = auxold(ifine,jfine,m)
                  endif
               enddo
               do m = 1,meqn
                  ql(ifine,m) = qold_coarse(ic,jc,m)
                  if (ifine .gt. -1) then
                     qr(ifine-1,m) = qold(ifine,jfine,m)
                  endif
               enddo
            enddo
         enddo

         ixy = 2
         call rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &         wave,s,amdq,apdq)

         do m = 1,meqn
            do ic = 1,mxc
               do ir = 1,iref_ratio
                  ifine = (ic-1)*iref_ratio + ir
                  if (ifine .lt. -1 .or. ifine .gt. mx+mbc) cycle
                  qadd_y(ifine,my+1,m) = (amdq(ifine,m) + apdq(ifine,m))
               enddo
            enddo
         enddo
      endif

      inqad = .false.

      end
