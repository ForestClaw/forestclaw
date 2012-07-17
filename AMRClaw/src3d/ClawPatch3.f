      subroutine ClawPatch3(maxm, meqn, maux, mbc, method, mthlim,
     &      mcapa, mwaves, mx, my, mz, qold, aux, dx, dy, dz, dt, cfl,
     &      work, mwork, qold_coarse, aux_coarse,
     &      qadd_x, qadd_y, qadd_z, auxtype_int, xlower, ylower, zlower,
     &      intersectsBoundary, level,mthbc,t,mxc, myc, mzc,
     &      fp,fm, gp, gm, hp, hm,
     &      fp_chombo, fm_chombo, gp_chombo,gm_chombo,
     &      hp_chombo, hm_chombo, fpc_chombo, fmc_chombo, gpc_chombo,
     &      gmc_chombo, hpc_chombo, hmc_chombo)

      implicit none
      external rpn3,rpt3, rptt3

      integer maxm,meqn,maux,mbc,mcapa,mwaves,mx,my, mz, mwork,
     &      maxmx, maxmy, maxmz, mxc, myc, mzc,level,
     &      intersectsBoundary(6),mthbc(6), mthbc_amr(6)
      double precision dx,dy,dz,dt,cfl, xlower, ylower, zlower,
     &      t, dtdx, dtdy, dtdz

      integer method(7), mthlim(mwaves)
      double precision work(mwork)

      double precision   qold(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &      1-mbc:mz+mbc,meqn)
      double precision   qold_coarse(0:mxc+1, 0:myc+1,
     &      0:mzc+1, meqn)

      double precision  aux(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &      1-mbc:mz+mbc, maux)

      double precision   aux_coarse(0:mxc+1, 0:myc+1,
     &      0:mzc+1, maux)

      integer auxtype_int(maux)
      integer maxmaux
      parameter(maxmaux = 100)
      character*10 auxtype(maxmaux)

c     Conservative numerical fluxes returned here.
      double precision fp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision fm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision gp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision gm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision hp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision hm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

      double precision fp_chombo(mx+1,my,  mz,   meqn)
      double precision fm_chombo(mx+1,my,  mz,   meqn)
      double precision gp_chombo(mx,  my+1,mz,   meqn)
      double precision gm_chombo(mx,  my+1,mz,   meqn)
      double precision hp_chombo(mx,  my,  mz+1, meqn)
      double precision hm_chombo(mx,  my,  mz+1, meqn)

      double precision fpc_chombo(mx+1,my,  mz,   meqn)
      double precision fmc_chombo(mx+1,my,  mz,   meqn)
      double precision gpc_chombo(mx,  my+1,mz,   meqn)
      double precision gmc_chombo(mx,  my+1,mz,   meqn)
      double precision hpc_chombo(mx,  my,  mz+1, meqn)
      double precision hmc_chombo(mx,  my,  mz+1, meqn)

      double precision qadd_x(mx+1,my,  mz,  meqn)
      double precision qadd_y(mx,  my+1,mz,  meqn)
      double precision qadd_z(mx,  my,  mz+1,meqn)

c     # Local variables
      integer i0faddm, i0faddp, i0gadd, i0hadd,
     &      i0q1d, i0dtdx1, i0dtdy1,i0dtdz1,
     &      i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl,
     &      i0auxr
      integer i,j,k,m, ic, jc, kc, refratio

c     # For debugging
      integer iunit, some_number
      logical debug

c     Needed by Riemann solvers.
      double precision dtcom, dxcom,dycom,dzcom, tcom
      integer icom, jcom, kcom
      common/comxyt/dtcom,dxcom,dycom,dzcom, tcom,icom,jcom, kcom

c -----------------------------------------------------------------
c     # Fix auxtype - convert from integers to characters
c     # This is a hack, but I didn't want to figure out how to
c     # convert std:vector<string> to arrays of character*10 needed by fortran
C     ...(DAC)

      if (maux .gt. maxmaux) then
         write(6,*) 'ClawPatch3 : maxmaux > maux;  Increase ',
     &         'size of maxmaux'
         stop
      endif

      do m = 1,maux
         if (auxtype_int(m) .eq. 1) then
            auxtype(m) = "xleft"
         elseif (auxtype_int(m) .eq. 2) then
            auxtype(m) = "yleft"
         elseif (auxtype_int(m) .eq. 3) then
            auxtype(m) = "zleft"
         elseif (auxtype_int(m) .eq. 4) then
            auxtype(m) = "center"
         elseif (auxtype_int(m) .eq. 5) then
            auxtype(m) = "capacity"
         endif
      enddo

c ------------------------------------------------------------------
c     # Set variables in common block, needed by step2.f and flux2.f
      call set_common_step3_flux3(method,mwaves,mthlim,mcapa)

c     This should be set to actual time, in case the user wants it
c     it for some reason in the Riemann solver.
      tcom = t
      dxcom = dx
      dycom = dy
      dzcom = dz

c     Not totally clear why these are set to zero... but it was done
c     this way in original AMRclaw code.
      method(1) = 0
      method(4) = 0

      maxmx = mx
      maxmy = my
      maxmz = mz

c -----------------------------------------------------------------------
c Some debugging information for iunit.
      debug = .false.
      if (debug) then
         iunit = 51
         open(iunit,file='clawpatch3.out')
         write(iunit,'(A,I5)') 'mx     = ', mx
         write(iunit,'(A,I5)') 'my     = ', my
         write(iunit,'(A,I5)') 'mz     = ', mz
         write(iunit,'(A,I5)') 'mxc    = ', mxc
         write(iunit,'(A,I5)') 'myc    = ', myc
         write(iunit,'(A,I5)') 'mzc    = ', mzc
         write(iunit,'(A,I5)') 'meqn   = ', meqn
         write(iunit,'(A,I5)') 'mwaves = ', mwaves
         write(iunit,'(A,I5)') 'mbc    = ', mbc
         write(iunit,'(A,I5)') 'maux   = ', maux
         write(iunit,'(A,I5)') 'mcapa  = ', mcapa
      endif

c -----------------------------------------------------------------------
c     # Set up boundary conditions
      do i = 1,6
         if (intersectsBoundary(i) .eq. 0) then
            mthbc_amr(i) = -1
         else
            mthbc_amr(i) = mthbc(i)
         endif
      enddo

      if (debug) then
c         write(iunit,'(A,4I4)') 'mthbc_tmp: ',(mthbc_amr(i),i=1,4)
      endif

      call bc3(maxmx, maxmy, maxmz, meqn, mbc, mx, my, mz,
     &      xlower, ylower, zlower,
     &      dx, dy, dz, qold, maux, aux, t, dt, mthbc_amr)

c ----------------------------------------------------------------------
c     # Before qold is modified, handle non-conservative fix-up at
c     # coarse/fine boundaries by modifying qadd_coarse.
      if (level .ne. 0) then
         i0wave  = 1
         i0s     = i0wave +  (maxm+2*mbc)*meqn*mwaves
         i0amdq  = i0s +     (maxm+2*mbc)*mwaves
         i0apdq  = i0amdq +  (maxm+2*mbc)*meqn
         i0ql    = i0apdq +  (maxm+2*mbc)*meqn
         i0qr    = i0ql +    (maxm+2*mbc)*meqn
         i0auxl  = i0qr +    (maxm+2*mbc)*meqn
         i0auxr  = i0auxl +  (maxm+2*mbc)*maux
         mused   = i0auxr +  (maxm+2*mbc)*maux - 1
c
         if (mused.gt.mwork) then
            write(6,*) 'Not enough work space in ClawPatch3'
            write(6,*) 'mused = ', mused, '   mwork =',mwork
            stop
         endif


         call qad3(maxm,meqn,mwaves,mbc,mx,my,mz,maux,qold,
     &         qold_coarse, aux, aux_coarse,
     &         intersectsBoundary, mxc,myc,mzc, dt,dx,dy,dz,
     &         work(i0ql),work(i0qr), work(i0auxl),work(i0auxr),
     &         work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &         qadd_x,qadd_y,qadd_z,auxtype)
      endif
c -----------------------------------------------------------------------
c     # take one step on the conservation law:


c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
      i0faddm = 1
      i0faddp = i0faddm +   (maxm+2*mbc)*meqn
      i0gadd  = i0faddp +   (maxm+2*mbc)*meqn
      i0hadd  = i0gadd  + 6*(maxm+2*mbc)*meqn
      i0q1d   = i0hadd  + 6*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d   +   (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1 +   (maxm+2*mbc)
      i0dtdz1 = i0dtdy1 +   (maxm+2*mbc)
      i0aux1 = i0dtdz1  +   (maxm+2*mbc)
      i0aux2 = i0aux1   +   (maxm+2*mbc)*maux*3
      i0aux3 = i0aux2   +   (maxm+2*mbc)*maux*3
c
c
      i0next = i0aux3 + (maxm+2*mbc)*maux*3    !# next free space
      mused  = i0next - 1                    !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step3)
c ----------------------------

      if (mused.gt.mwork) then
c        # This shouldn't happen due to checks in claw3
         write(6,*) 'Not enough work space in clawpatch3'
         write(6,*) 'mused = ', mused, '   mwork =',mwork
         stop
      endif

      call b4step3(mx,my,mz, mbc,mx,my,mz, meqn,qold,
     &      xlower,ylower,zlower, dx,dy,dz, t,dt,maux,aux)

c
c
c     # take one step on the conservation law:
c
      call step3(maxm,maxmx,maxmy,maxmz, meqn,maux, mbc,
     &      mx,my, mz, qold,aux,dx,dy,dz, dt,
     &      cfl,fm,fp,gm,gp,hm, hp,
     &      work(i0faddm),work(i0faddp),
     &      work(i0gadd),work(i0hadd),
     &      work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &      work(i0dtdz1),
     &      work(i0aux1),work(i0aux2),work(i0aux3),
     &      work(i0next),mwork1,rpn3,rpt3,rptt3)


c
c
c ---------------------------------------------------------------

c     # update q
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
      do m = 1,meqn
         do i = 1,mx
            do j = 1,my
               do k = 1,mz
                  if (mcapa.eq.0) then
c                    # no capa array.  Standard flux differencing:
                     qold(i,j,k,m) = qold(i,j,k,m)
     &                     - dtdx * (fm(i+1,j,k,m) - fp(i,j,k,m))
     &                     - dtdy * (gm(i,j+1,k,m) - gp(i,j,k,m))
     &                     - dtdz * (hm(i,j,k+1,m) - hp(i,j,k,m))
                  else
                     qold(i,j,k,m) = qold(i,j,k,m)
     &                     - (dtdx * (fm(i+1,j,k,m) - fp(i,j,k,m))
     &                     +  dtdy * (gm(i,j+1,k,m) - gp(i,j,k,m))
     &                     +  dtdz * (hm(i,j,k+1,m) - hp(i,j,k,m)))
     &                     / aux(i,j,k,mcapa)
                  endif
               enddo
            enddo
         enddo
      enddo

      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
c        # In AMRClaw, a function called 'src1d' is called from
c        # from qad.  But we don't do that here.
         call src3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &         xlower,ylower,zlower,dx,dy,dz,
     &         qold,maux,aux,t,dt)
      endif


c ---------------------------------------------------------------
c
c     # Assign flux registers fp_chombo, etc.  The reason we do this, rather
c     # then pass fp_chombo, etc directly into step2, is that fp and
c     # fp_chombo (fm and fm_chombo, etc.) are dimensioned differently.
c

       refratio = mx/mxc
       if (mcapa .gt. 0) then
          do i = 1,mx+1
             ic = (i-1)/refratio + 1
             do j = 1,my
                jc = (j-1)/refratio + 1
                do k = 1,mz
                   kc = (k-1)/refratio + 1
                   do m = 1,meqn
                      fp_chombo(i,j,k,m) =  fp(i,j,k,m) /
     &                      aux(i,j,k,mcapa)

                      fm_chombo(i,j,k,m) = fm(i,j,k,m) /
     &                      aux(i-1,j,k,mcapa)

                      fpc_chombo(i,j,k,m) = fp(i,j,k,m) /
     &                      aux_coarse(ic,jc,kc,mcapa)

                      fmc_chombo(i,j,k,m) = fm(i,j,k,m) /
     &                      aux_coarse(ic-1,jc,kc,mcapa)

                      if (i .eq. 1) then
                         qadd_x(i,j,k,m) = qadd_x(i,j,k,m) /
     &                         aux_coarse(ic-1,jc,kc,mcapa)
                      elseif (i .eq. mx+1) then
                         qadd_x(i,j,k,m) = qadd_x(i,j,k,m) /
     &                         aux_coarse(ic,jc,kc,mcapa)
                      endif
                   enddo
                enddo
             enddo
          enddo

          do i = 1,mx
             ic = (i-1)/refratio + 1
             do j = 1,my+1
                jc = (j-1)/refratio + 1
                do k = 1,mz
                   kc = (k-1)/refratio + 1
                   do m = 1,meqn
                      gp_chombo(i,j,k,m) = gp(i,j,k,m)
     &                      /aux(i,j,k,mcapa)
                      gm_chombo(i,j,k,m) = gm(i,j,k,m)/
     &                      aux(i,j-1,k,mcapa)

                      gpc_chombo(i,j,k,m) = gp(i,j,k,m) /
     &                      aux_coarse(ic,jc,kc,mcapa)

                      gmc_chombo(i,j,k, m) = gm(i,j,k,m)/
     &                      aux_coarse(ic,jc-1,kc, mcapa)

                      if (j .eq. 1) then
                         qadd_y(i,j,k,m) = qadd_y(i,j,k,m) /
     &                         aux_coarse(ic,jc-1,kc,mcapa)
                      elseif (j .eq. my+1) then
                         qadd_y(i,j,k,m) = qadd_y(i,j,k,m) /
     &                         aux_coarse(ic,jc,kc,mcapa)
                      endif
                   enddo
                enddo
             enddo
          enddo

          do i = 1,mx
             ic = (i-1)/refratio + 1
             do j = 1,my
                jc = (j-1)/refratio + 1
                do k = 1,mz+1
                   kc = (k-1)/refratio + 1
                   do m = 1,meqn
                      hp_chombo(i,j,k,m) = hp(i,j,k,m)
     &                      /aux(i,j,k,mcapa)
                      hm_chombo(i,j,k,m) = hm(i,j,k,m)/
     &                      aux(i,j,k-1,mcapa)

                      hpc_chombo(i,j,k,m) = hp(i,j,k,m) /
     &                      aux_coarse(ic,jc,kc,mcapa)

                      hmc_chombo(i,j,k, m) = hm(i,j,k,m)/
     &                      aux_coarse(ic,jc,kc-1, mcapa)

                      if (k .eq. 1) then
                         qadd_z(i,j,k,m) = qadd_z(i,j,k,m) /
     &                         aux_coarse(ic,jc,kc-1,mcapa)
                      elseif (k .eq. mz+1) then
                         qadd_z(i,j,k,m) = qadd_z(i,j,k,m) /
     &                         aux_coarse(ic,jc,kc,mcapa)

                      endif
                   enddo
                enddo
             enddo
          enddo
       else
c         # no capacity function
          do i = 1,mx+1
             do j = 1,my
                do k = 1,mz
                   do m = 1,meqn
                      fp_chombo(i,j,k,m) = fp(i,j,k,m)
                      fm_chombo(i,j,k,m) = fm(i,j,k,m)
                      fpc_chombo(i,j,k,m) = fp(i,j,k,m)
                      fmc_chombo(i,j,k,m) = fm(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo

          do i = 1,mx
             do j = 1,my+1
                do k = 1,mz
                   do m = 1,meqn
                      gp_chombo(i,j,k,m) = gp(i,j,k,m)
                      gm_chombo(i,j,k,m) = gm(i,j,k,m)
                      gpc_chombo(i,j,k,m) = gp(i,j,k,m)
                      gmc_chombo(i,j,k,m) = gm(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo

          do i = 1,mx
             do j = 1,my
                do k = 1,mz+1
                   do m = 1,meqn
                      hp_chombo(i,j,k,m) = hp(i,j,k,m)
                      hm_chombo(i,j,k,m) = hm(i,j,k,m)
                      hpc_chombo(i,j,k,m) = hp(i,j,k,m)
                      hmc_chombo(i,j,k,m) = hm(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       endif

c      ---------------------------------------------------------

      if (debug) then
         write(iunit,*) 'qold...'
         do i = 1,mx
            do j = 1,my
               write(iunit,*) 'i = ', i, ', j = ', j
               do k = 1,mz
                  write(iunit,'(3I5,4E24.16)') i,j,k,qold(i,j,k,1),
     &                  fm(i,j,k,1), fp(i,j,k,1)
               enddo
               write(iunit,*) ' '
            enddo
            write(iunit,*) ' '
         enddo
      endif

      if (debug) then
         close(iunit)
      endif

      if (debug) then
         write(6,*) 'In Clawpatch3; input a number to continue'
         read(5,*) some_number
      endif

      end


c ----------------------------------------------------------------------
c This sets up the common blocks needed by step3.f and flux3.f. It
c replaces the include "call.i" in these routines. No other changes
c to these routines is required.
c This routine is identical to set_common_step2_flux2, in ClawPatch2.f
c ----------------------------------------------------------------------

      subroutine set_common_step3_flux3(a_method,a_mwaves,
     &      a_mthlim,a_mcapa)
      implicit none

      include "claw.i"

      integer a_mwaves, a_mcapa
      integer a_method(7), a_mthlim(a_mwaves)
      integer i

      mcapa = a_mcapa
      mwaves = a_mwaves

      do i = 1,7
         method(i) = a_method(i)
      enddo

      if (mwaves .gt. max_mwaves) then
         write(6,*) 'set_common3 : mwaves > max_mwaves (in claw.i)'
         stop
      endif

      do i = 1,mwaves
         mthlim(i) = a_mthlim(i)
      enddo

      end


c     # This is called from ClawPatch.cpp, in ClawPatch::setAuxArray().
c     # this same subroutine is also in src2d
      subroutine set_common_levels(a_maxlevel,a_level,a_refratios)
      implicit none

c     # Inputs
      integer a_level, a_maxlevel, i
      integer a_refratios(a_maxlevel)

c     # set common block  that can be seen by setaux.f, for example.
      integer max_maxlevel
      parameter(max_maxlevel = 100)
      integer level, refratios(max_maxlevel), maxlevel
      common /comlevel/ maxlevel, level, refratios

      if (a_maxlevel .gt. max_maxlevel) then
         write(6,'(A,A)') 'set_common_levels : ',
     &         'maxlevel >= max_maxlevel '
         stop
      endif

      maxlevel = a_maxlevel  !! maximum grid level; maxlevel = 0 --> no refinement

      level = a_level  !! Current level

      do i = 1,maxlevel
         refratios(i) = a_refratios(i)
      enddo

      end
