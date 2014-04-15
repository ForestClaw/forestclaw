      subroutine clawpatch2(maxm, meqn, maux, mbc, method, mthlim,
     &      mcapa, mwaves, mx, my, qold, aux, dx, dy, dt, cfl,
     &      xlower,ylower,level,t, fp,fm, gp, gm)

      implicit none

      external rpn2,rpt2

      integer maxm,meqn,maux,mbc,mcapa,mwaves,mx,my, mwork,
     &      maxmx, maxmy, level

      double precision dx,dy,dt,cfl, xlower, ylower, t,
     &      dtdx, dtdy

      integer method(7), mthlim(mwaves)
      integer maxmaux
      parameter(maxmaux = 100)
c      double precision work(mwork)

      double precision qold(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

c     Conservative numerical fluxes returned here.
      double precision fp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision fm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision gp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision gm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)


cc     # Local variables
c      integer i0faddm, i0faddp, i0gaddm, i0gaddp,
c     &      i0q1d, i0dtdx1, i0dtdy1,
c     &      i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
c
c      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl,
c     &      i0auxr
      integer i,j,m, refratio, ic, jc

      integer iunit
      logical debug

c     Needed by Riemann solvers.
      double precision dtcom, dxcom,dycom,tcom
      integer icom, jcom
      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      data debug /.false./

c -----------------------------------------------------------------

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
c      i0faddm = 1
c      i0faddp = i0faddm +   (maxm+2*mbc)*meqn
c      i0gaddm = i0faddp +   (maxm+2*mbc)*meqn
c      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
c      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
c      i0dtdx1 = i0q1d   +   (maxm+2*mbc)*meqn
c      i0dtdy1 = i0dtdx1 +   (maxm+2*mbc)
c      i0aux1  = i0dtdy1 +   (maxm+2*mbc)
c      i0aux2  = i0aux1  +   (maxm+2*mbc)*maux
c      i0aux3  = i0aux2  +   (maxm+2*mbc)*maux
cc
cc
c      i0next  = i0aux3  + (maxm+2*mbc)*maux    !# next free space
c      mused   = i0next - 1                    !# space already used
c      mwork1  = mwork - mused           !# remaining space (passed to step2)

c      if (mused.gt.mwork) then
cc        # This shouldn't happen due to checks in claw2
c         write(6,*) 'Not enough work space in clawpatch2'
c         write(6,*) 'mused = ', mused, '   mwork =',mwork
c         stop
c      endif

c     # ----------------------------------------------
c     # Take one step on the conservation law:
c     # ----------------------------------------------

      call  step2(maxm,meqn,maux,mbc,mx,my,qold,aux,dx,dy,dt,
     &      cfl,fm,fp,gm,gp,rpn2,rpt2)

c     # ----------------------------------------------
c     # Update q
c     # ----------------------------------------------
      dtdx = dt/dx
      dtdy = dt/dy
      do m = 1,meqn
         do i = 1,mx
            do j = 1,my
               if (mcapa.eq.0) then
c                 # no capa array.  Standard flux differencing:
                  qold(m,i,j) = qold(m,i,j)
     &                  - dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &                  - dtdy * (gm(m,i,j+1) - gp(m,i,j))
               else
c                 # with capa array.
                  qold(m,i,j) = qold(m,i,j)
     &                  -(dtdx*(fm(m,i+1,j) - fp(m,i,j))
     &                  + dtdy*(gm(m,i,j+1) - gp(m,i,j)))/aux(mcapa,i,j)
               endif
            enddo
         enddo
      enddo


c      if (method(5).eq.1) then
cc        # with source term:   use Godunov splitting
cc        # In AMRClaw, a function called 'src1d' is called from
cc        # from qad.  But we don't do that here.
cc
cc        # Not clear where this should fit in now.
c         call src2(maxmx,maxmy,meqn,mbc,mx,my,
c     &         xlower,ylower,dx,dy,
c     &         qold,maux,aux,t,dt)
c      endif
c
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
