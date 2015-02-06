      subroutine clawpatch2(maxm, meqn, maux, mbc, method, mthlim,
     &      mcapa, mwaves, mx, my, qold, aux, dx, dy, dt, cfl,
     &      work, mwork,xlower,ylower,level,t, fp,fm, gp, gm,
     &      rpn2, rpt2)

      implicit none

      external rpn2,rpt2

      integer maxm,meqn,maux,mbc,mcapa,mwaves,mx,my, mwork,
     &      maxmx, maxmy, level
      double precision dx,dy,dt,cfl, xlower, ylower, t,
     &      dtdx, dtdy

      integer method(7), mthlim(mwaves)
      integer maxmaux
      parameter(maxmaux = 100)
      double precision work(mwork)

      double precision qold(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

c     Conservative numerical fluxes returned here.
      double precision fp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision fm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)


c     # Local variables
      integer i0faddm, i0faddp, i0gaddm, i0gaddp,
     &      i0q1d, i0dtdx1, i0dtdy1,
     &      i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl,
     &      i0auxr
      integer i,j,m, refratio, ic, jc

      integer iunit
      logical debug

c     Needed by Riemann solvers.
      double precision dtcom, dxcom,dycom,tcom
      integer icom, jcom
      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      data debug /.false./

c -----------------------------------------------------------------
c     # Fix auxtype - convert from integers to characters
c     # This is a hack, but I didn't want to figure out how to
c     # convert std:vector<string> to arrays of character*10 needed by fortran
C     ...(DAC)

      if (maux .gt. maxmaux) then
         write(6,*) 'ClawPatch2 : maxmaux > maux;  Increase ',
     &         'size of maxmaux'
         stop
      endif

c ------------------------------------------------------------------
c     # Set variables in common block, needed by step2.f and flux2.f
      call set_common_step2_flux2(method,mwaves,mthlim,mcapa)

c     This should be set to actual time, in case the user wants it
c     it for some reason in the Riemann solver.
      tcom = t
      dxcom = dx
      dycom = dy

c     Not totally clear why these are set to zero... but it was done
c     this way in original AMRclaw code.
      method(1) = 0
      method(4) = 0

      maxmx = mx
      maxmy = my

c -----------------------------------------------------------------------
c Some debugging information for iunit.
      debug = .false.
      if (debug) then
         iunit = 51
         open(iunit,file='clawpatch2.out')
         write(iunit,'(A,I5)') 'mx     = ', mx
         write(iunit,'(A,I5)') 'my     = ', my
         write(iunit,'(A,I5)') 'meqn   = ', meqn
         write(iunit,'(A,I5)') 'mwaves = ', mwaves
         write(iunit,'(A,I5)') 'mbc    = ', mbc
         write(iunit,'(A,I5)') 'maux   = ', maux
         write(iunit,'(A,I5)') 'mcapa  = ', mcapa
      endif

c -----------------------------------------------------------------------
c     # take one step on the conservation law:

c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
      i0faddm = 1
      i0faddp = i0faddm +   (maxm+2*mbc)*meqn
      i0gaddm = i0faddp +   (maxm+2*mbc)*meqn
      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d   +   (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1 +   (maxm+2*mbc)
      i0aux1  = i0dtdy1 +   (maxm+2*mbc)
      i0aux2  = i0aux1  +   (maxm+2*mbc)*maux
      i0aux3  = i0aux2  +   (maxm+2*mbc)*maux
c
c
      i0next  = i0aux3  + (maxm+2*mbc)*maux    !# next free space
      mused   = i0next - 1                    !# space already used
      mwork1  = mwork - mused           !# remaining space (passed to step2)

      if (mused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) 'Not enough work space in clawpatch2'
         write(6,*) 'mused = ', mused, '   mwork =',mwork
         stop
      endif

c     # This is now called from C code
c      call b4step2(mx,my,mbc,mx,my,meqn,qold, xlower,ylower,dx,dy,t,
c     &      dt,maux,aux)

c
c     # take one step on the conservation law:
c
      call step2(maxm,maxmx,maxmy,meqn,maux, mbc,
     &      mx,my, qold,aux,dx,dy,dt,
     &      cfl,fm,fp,gm,gp,
     &      work(i0faddm),work(i0faddp),
     &      work(i0gaddm),work(i0gaddp),
     &      work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &      work(i0aux1),work(i0aux2),work(i0aux3),
     &      work(i0next),mwork1,rpn2,rpt2)

c
c
c ---------------------------------------------------------------
c
c     # update q
      dtdx = dt/dx
      dtdy = dt/dy
      do m = 1,meqn
         do i = 1,mx
            do j = 1,my
               if (mcapa.eq.0) then
c                 # no capa array.  Standard flux differencing:
                  qold(i,j,m) = qold(i,j,m)
     &                  - dtdx * (fm(i+1,j,m) - fp(i,j,m))
     &                  - dtdy * (gm(i,j+1,m) - gp(i,j,m))
               else
c                 # with capa array.
                  qold(i,j,m) = qold(i,j,m)
     &                  -(dtdx*(fm(i+1,j,m) - fp(i,j,m))
     &                  + dtdy*(gm(i,j+1,m) - gp(i,j,m)))/aux(i,j,mcapa)
               endif
            enddo
         enddo
      enddo

      end


c     # This sets variables for step2 and flux2 that are not passed in
c     # through the argument list.  This is for compatibility with amrclaw,
c     # and not because of some grand software design issue...
      subroutine set_common_step2_flux2(a_method,a_mwaves,
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
         write(6,*) 'set_common2 : mwaves > max_mwaves (in claw.i)'
         stop
      endif

      do i = 1,mwaves
         mthlim(i) = a_mthlim(i)
      enddo

      end

cc     # This is called from ClawPatch.cpp, in ClawPatch::setAuxArray().
cc     # this same subroutine is also in src2d
c      subroutine set_common_levels(a_maxlevel,a_level,a_refratio)
c      implicit none
c
cc     # Inputs
c      integer a_level, a_maxlevel, a_refratio
c
cc     # set common block  that can be seen by setaux.f, for example.
c      integer com_level, com_maxlevel, com_refratio
c      common /comlevel/ com_maxlevel, com_level, com_refratio
c
c      com_maxlevel = a_maxlevel  !! maximum grid level; maxlevel = 0 --> no refinement
c      com_level = a_level  !! Current level
c      com_refratio = a_refratio
c
c      end
