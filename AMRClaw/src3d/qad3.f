      subroutine qad3(maxm,meqn,mwaves,mbc,mx,my,mz, maux,
     &      qold,qold_coarse, auxold, auxold_coarse,
     &      intersectsBoundary,
     &      mxc,myc,mzc, dt,dx,dy, dz,ql,qr,auxl,auxr, wave,s,
     &      amdq,apdq, qadd_x,qadd_y,qadd_z,auxtype)


      implicit none
c     # Input parameters
      integer maxm,meqn,mwaves,mbc,mx,my,mz,mxc,myc,mzc,maux
      integer intersectsBoundary(6)
      double precision dt, dx, dy,dz

c     # note that mbc is only relevant for the fine grid;  the coarse
c     # grid only has a single layer of ghost cells.
      double precision qold(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,
     &      meqn)
      double precision qold_coarse(0:mxc+1,0:myc+1,0:mzc+1,meqn)

      double precision auxold(1-mbc:mx+mbc,1-mbc:my+mbc,
     &      1-mbc:mz+mbc,maux)
      double precision auxold_coarse(0:mxc+1,0:myc+1,0:mzc+1,maux)

      double precision qadd_x(mx+1,my,  mz,  meqn)
      double precision qadd_y(mx,  my+1,mz,  meqn)
      double precision qadd_z(mx,  my  ,mz+1,meqn)

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)

      double precision auxl(1-mbc:maxm+mbc,maux)
      double precision auxr(1-mbc:maxm+mbc,maux)

      character*10 auxtype(maux)

      integer ic,jc,kc, iff,jf,kf, m, iref_ratio, i, ixyz
      logical ib(6)
      integer ir,jr,kr

c     For debugging
      integer iunit
      logical debug, inqad

      common /comqad/ inqad

      inqad = .true.

      debug = .false.
      if (debug) then
         iunit = 55
         open(iunit, file='qad.out')
         write(iunit,'(A,i5)') 'meqn = ', meqn
         write(iunit,'(A,i5)') 'mwaves = ', mwaves
         write(iunit,'(A,I5)') 'mx = ', mx
         write(iunit,'(A,i5)') 'my = ', my
         write(iunit,'(A,i5)') 'mz = ', mz
         write(iunit,'(A,I5)') 'mxc = ', mxc
         write(iunit,'(A,i5)') 'myc = ', myc
         write(iunit,'(A,i5)') 'mzc = ', mzc
         write(iunit,'(A,I5)') 'maxm = ', maxm
         write(iunit,'(A,I5)') 'refratio = ', mx/mxc
         write(iunit,'(A,6I5)') 'bdry info : ',
     &         (intersectsBoundary(m),m=1,6)
      endif

      iref_ratio = mx/mxc

c     # Now check boundaries to get something reasonable if the
c     # coarse grid intersects a physical boundary.
c     # Note that it doesn't really matter what we put here, but
c     # we want something that makes sense for the problem so that
c     # when we pass everything to the Riemann solver, we don't
c     # crash with non-sensical values.

      do i = 1,6
         ib(i) = intersectsBoundary(i) .eq. 1
      enddo

      if (ib(1)) then
         do  jc = 0,myc+1
            do kc = 0,mzc+1
               do m = 1,meqn
                  qold_coarse(0,jc,kc,m) = qold_coarse(1,jc,kc,m)
               enddo
               do m = 1,maux
                  auxold_coarse(0,jc,kc,m) = auxold_coarse(1,jc,kc,m)
               enddo
            enddo
         enddo
      endif

      if (ib(2)) then
         do  jc = 0,myc+1
            do kc = 0,mzc+1
               do m = 1,meqn
                  qold_coarse(mxc+1,jc,kc, m) = qold_coarse(mxc,jc,kc,m)
               enddo
               do m = 1,maux
                  auxold_coarse(mxc+1,jc,kc, m) =
     &                  auxold_coarse(mxc,jc,kc,m)
               enddo
            enddo
         enddo
      endif


      if (ib(3)) then
         do  ic = 0,mxc+1
            do kc = 0,mzc+1
               do m = 1,meqn
                  qold_coarse(ic,0,kc,m) = qold_coarse(ic,1,kc,m)
               enddo
               do m = 1,maux
                  auxold_coarse(ic,0,kc,m) = auxold_coarse(ic,1,kc,m)
               enddo
            enddo
         enddo
      endif

      if (ib(4)) then
         do  ic = 0,mxc+1
            do kc = 0,mzc+1
               do m = 1,meqn
                  qold_coarse(ic,myc+1,kc,m) = qold_coarse(ic,myc,kc,m)
               enddo
               do m = 1,maux
                  auxold_coarse(ic,myc+1,kc,m) =
     &                  auxold_coarse(ic,myc,kc,m)
               enddo
            enddo
         enddo
      endif


      if (ib(5)) then
         do ic = 0,mxc+1
            do jc = 0,myc+1
               do m = 1,meqn
                  qold_coarse(ic,jc,0,m) = qold_coarse(ic,jc,1,m)
               enddo
               do m = 1,maux
                  auxold_coarse(ic,jc,0,m) = auxold_coarse(ic,jc,1,m)
               enddo
            enddo
         enddo
      endif

      if (ib(6)) then
         do ic = 0,mxc+1
            do jc = 0,myc+1
               do m = 1,meqn
                  qold_coarse(ic,jc,mzc+1,m) = qold_coarse(ic,jc,mzc,m)
               enddo
               do m = 1,maux
                  auxold_coarse(ic,jc,mzc+1,m) =
     &                  auxold_coarse(ic,jc,mzc,m)
               enddo
            enddo
         enddo
      endif

      if (iref_ratio .eq. 1) then
c        # Should be equvilient to
c        #           (ib(1) .and. ib(2) .and. ib(3) .and. ib(4))
c        #
c        # We are on the coarsest grid,and there is really nothing we
c        # should be doing here, except perhaps zeroing out qadd_coarse
c        # which we did above.
         if (debug) then
            write(iunit,*) 'Returning from qad, since ref_ratio == 1'
            close(iunit)
         endif
         return
      endif

c ----------------------------------------------------------------
c     # Side 1
c     # Left edge (xleft)
c ----------------------------------------------------------------

      if (.not. ib(1)) then
         if (debug) then
            write(iunit,*) 'Solving for side 1...'
            write(6,*) 'Solving for side 1...'
         endif
         ic = 0
         iff = 0
c        # Loop over coarse grid z values
         do kc = 1,mzc
            do kr = 1,iref_ratio
               kf = iref_ratio*(kc - 1) + kr
               do jc = 0,myc+1
                  do jr = 1,iref_ratio
                     jf = iref_ratio*(jc - 1) + jr
                     do m = 1,maux
                        if (auxtype(m) .eq. "xleft") then
                           auxl(jf,m) = auxold(iff+1,jf,kf,m)
                        else
                           auxl(jf,m) = auxold(iff,jf,kf,m)
                        endif
                        if (jf .gt. -1) then
                           auxr(jf-1,m) = auxold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(jf,m) = qold(iff,jf,kf,m)
                        if (jf .gt. -1) then
                           qr(jf-1,m) = qold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                  enddo
               enddo                    !! end meqn loop

c              # Doing a y sweep
               ixyz = 1
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,my,ql,qr,auxl,auxr,
     &               maux, wave,s,amdq,apdq)

               do m = 1,meqn
                  do jc = 1,myc
                     do jr = 1,iref_ratio
                        jf = iref_ratio*(jc-1) + jr
                        qadd_x(1,jf,kf,m) = (amdq(jf,m) + apdq(jf,m))
                     enddo
                  enddo
               enddo
            enddo                       !! end kr loop over fine grid z values
         enddo                          !! end loop kc over coarse grid z values
         if (debug) then
            write(6,*) 'Done with side 1...'
         endif
      endif


c ----------------------------------------------------------------
c     # Side 2
c     # left edge (xright)
c ----------------------------------------------------------------
      if (.not. ib(2)) then
         if (debug) then
            write(iunit,*) 'Solving for side 2...'
            write(6,*) 'Solving for side 2...'
         endif
         ic = mxc + 1
         iff = mx + 1
         do kc = 1,mzc
            do kr = 1,iref_ratio
               kf = iref_ratio*(kc - 1) + kr
c              # Do not loop 1-mbc, myc+mbc - the coarse grid only has
c              # a single layter of ghost cells.
c              # See note for side 1 as to why we loop over ghost cells
C              # here...
               do jc = 0,myc+1
c                 # Loop over number of cells in fine grid per coarse grid
C                 # cell
                  do jr = 1,iref_ratio
                     jf = iref_ratio*(jc - 1) + jr
                     do m = 1,maux
                        auxl(jf,m) = auxold_coarse(ic,jc,kc,m)
                        if (jf .gt. -1) then
                           auxr(jf-1,m) = auxold(iff,jf,kf,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(jf,m) = qold_coarse(ic,jc,kc,m)
                        if (jf .gt. -1) then
                           qr(jf-1,m) = qold(iff,jf,kf,m)
                        endif
                     enddo
                  enddo
               enddo

c              # Doing an x-sweep
               ixyz = 1
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,my,ql,qr,auxl,auxr,
     &               maux,wave,s,amdq,apdq)

               do m = 1,meqn
                  do jc = 1,myc
                     do jr = 1,iref_ratio
                        jf = iref_ratio*(jc-1) + jr
                        qadd_x(mx+1,jf,kf,m) = (amdq(jf,m) + apdq(jf,m))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         if (debug) then
            write(6,*) 'Done with side 2...'
         endif
      endif
c     ----------------------------------------------------------------
c     # Side 3
c     # Lower edge (yleft)
c ----------------------------------------------------------------
      if (.not. ib(3)) then
         if (debug) then
            write(6,*) 'Solving for side 3...'
            write(iunit,*) 'Solving for side 3...'
         endif
         jc = 0
         jf = 0
         do kc = 1,mzc
            do kr = 1,iref_ratio
               kf = iref_ratio*(kc-1) + kr
c              # Do not loop 1-mbc, mxc+mbc - the coarse grid only has
c              # a single layter of ghost cells.
c              # See note for side 1 as to why we loop over ghost cells
C              # here...
               do ic = 0,mxc+1
c                 # These values aren't really right, but they won't get
C                 # used.
                  do ir = 1,iref_ratio
                     iff = iref_ratio*(ic-1) + ir
                     do m = 1,maux
                        if (auxtype(m) .eq. "yleft") then
                           auxl(iff,m) = auxold(iff,jf+1,kf,m)
                        else
                           auxl(iff,m) = auxold(iff,jf,kf,m)
                        endif
                        if (iff .gt. -1) then
                           auxr(iff-1,m) = auxold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(iff,m) = qold(iff,jf,kf,m)
                        if (iff .gt. -1) then
                           qr(iff-1,m) = qold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                  enddo
               enddo

               ixyz = 2
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &               maux, wave,s,amdq,apdq)

               do m = 1,meqn
                  do ic = 1,mxc
                     do ir = 1,iref_ratio
                        iff = iref_ratio*(ic-1) + ir
                        qadd_y(iff,1,kf,m) = (amdq(iff,m) + apdq(iff,m))
                     enddo
                  enddo
               enddo
            enddo                       !! end kr loop on fine z direction
         enddo                          !! end kc loop on coarse z direction
         if (debug) then
            write(6,*) 'Done with side 3...'
         endif

      endif



c ----------------------------------------------------------------
c     # Side 4
c     # Upper edge (yright)
c ----------------------------------------------------------------
      if (.not. ib(4)) then
         if (debug) then
            write(iunit,*) 'Solving for side 4...'
            write(6,*) 'Solving for side 4'
         endif
         jc = myc+1
         jf = my+1
         do kc = 1,mzc
            do kr = 1,iref_ratio
               kf = iref_ratio*(kc - 1) + kr
               do ic = 0,mxc+1
c                 # Do not loop 1-mbc, mxc+mbc - the coarse grid only has
c                 # a single layter of ghost cells.
c                 # See note for side 1 as to why we loop over ghost cells
C                 here...
                  do ir = 1,iref_ratio
                     iff = iref_ratio*(ic-1) + ir
                     do m = 1,maux
                        auxl(iff,m) = auxold_coarse(ic,jc,kc,m)
                        if (iff .gt. -1) then
                           auxr(iff-1,m) = auxold(iff,jf,kf,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(iff,m) = qold_coarse(ic,jc,kc,m)
                        if (iff .gt. -1) then
                           qr(iff-1,m) = qold(iff,jf,kf,m)
                        endif
                     enddo
                  enddo
               enddo

               ixyz = 2
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &               maux,wave,s,amdq,apdq)

               do m = 1,meqn
                  do ic = 1,mxc
                     do ir = 1,iref_ratio
                        iff = iref_ratio*(ic-1) + ir
                        qadd_y(iff,my+1,kf,m) = (amdq(iff,m) +
     &                        apdq(iff,m))
                     enddo
                  enddo
               enddo
            enddo                       !! end kr loop over fine z grid
         enddo                          !! end kc loop over coarse z grid
         if (debug) then
            write(6,*) 'done with side 4'
         endif
      endif

c ----------------------------------------------------------------
c     # Side 5
c     # Lower edge (zleft)
c ----------------------------------------------------------------
      if (.not. ib(5)) then
         if (debug) then
            write(6,*) 'Solving for side 5...'
            write(iunit,*) 'Solving for side 5...'
         endif
         kc = 0
         kf = 0
         do jc = 1,myc
            do jr = 1,iref_ratio
               jf = iref_ratio*(jc - 1) + jr
               do ic = 0,mxc+1
c                 # Do not loop 1-mbc, mxc+mbc - the coarse grid only has
c                 # a single layter of ghost cells.
c                 # See note for side 1 as to why we loop over ghost cells
C                 here...
                  do ir = 1,iref_ratio
                     iff = iref_ratio*(ic-1) + ir
                     do m = 1,maux
                        if (auxtype(m) .eq. "zleft") then
                           auxl(iff,m) = auxold(iff,jf,kf+1,m)
                        else
                           auxl(iff,m) = auxold(iff,jf,kf,m)
                        endif
                        if (iff .gt. -1) then
                           auxr(iff-1,m) = auxold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(iff,m) = qold(iff,jf,kf,m)
                        if (iff .gt. -1) then
                           qr(iff-1,m) = qold_coarse(ic,jc,kc,m)
                        endif
                     enddo
                  enddo
               enddo

               ixyz = 3
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &               maux,wave,s,amdq,apdq)

               do m = 1,meqn
                  do ic = 1,mxc
                     do ir = 1,iref_ratio
                        iff = iref_ratio*(ic-1) + ir
                        qadd_z(iff,jf,1,m) = (amdq(iff,m) +
     &                        apdq(iff,m))
                      enddo
                  enddo
               enddo
            enddo                       !! end kr loop over fine z grid
         enddo                          !! end kc loop over coarse z grid
         if (debug) then
            write(6,*) 'Done with side 5...'
         endif
      endif

c ----------------------------------------------------------------
c     # Side 6
c     # Upper edge (zupper)
c ----------------------------------------------------------------
      if (.not. ib(6)) then
         if (debug) then
            write(6,*) 'Solving for side 6...'
            write(iunit,*) 'Solving for side 6...'
         endif
         kc = mzc+1
         kf = mz+1
         do jc = 1,myc
            do jr = 1,iref_ratio
               jf = iref_ratio*(jc - 1) + jr
               do ic = 0,mxc+1
c                 # Avoid sending bogus values to Riemann solver, if possible.
c                 # These values aren't really correct, but they won't get
C                 used.
                  do ir = 1,iref_ratio
                     iff = iref_ratio*(ic-1) + ir
                     do m = 1,maux
                        auxl(iff,m) = auxold_coarse(ic,jc,kc,m)
                        if (iff .gt. -1) then
                           auxr(iff-1,m) = auxold(iff,jf,kf,m)
                        endif
                     enddo
                     do m = 1,meqn
                        ql(iff,m) = qold_coarse(ic,jc,kc,m)
                        if (iff .gt. -1) then
                           qr(iff-1,m) = qold(iff,jf,kf,m)
                        endif
                     enddo
                  enddo
               enddo

               ixyz = 3
               call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &               maux, wave,s,amdq,apdq)

               do m = 1,meqn
                  do ic = 1,mxc
                     do ir = 1,iref_ratio
                        iff = iref_ratio*(ic-1) + ir
                        qadd_z(iff,jf,mz+1,m) = (amdq(iff,m) +
     &                        apdq(iff,m))
                     enddo
                  enddo
               enddo
            enddo                       !! end kr loop over fine z grid
         enddo                          !! end kc loop over coarse z grid
      endif


c     # in amrclaw, a function src1d is also called from here, but I haven't
c     # put that in yet.

      inqad = .false.


      end
