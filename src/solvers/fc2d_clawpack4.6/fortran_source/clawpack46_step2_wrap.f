      subroutine clawpack46_step2_wrap(maxm, meqn, maux, mbc,
     &      method, mthlim, mcapa, mwaves, mx, my, qold, aux,
     &      dx, dy, dt,cfl, work, mwork,xlower,ylower,level,
     &      t, fp,fm, gp, gm, rpn2, rpt2, rpn2fw, rpt2fw, flux2,
     &      block_corner_count,ierror, use_fwaves)

      implicit none

      external rpn2,rpt2, rpn2fw, rpt2fw, flux2

      integer maxm,meqn,maux,mbc,mcapa,mwaves,mx,my, mwork
      integer maxmx, maxmy, level, ierror, use_fwaves
      integer method(7), mthlim(mwaves)
      integer block_corner_count(0:3)

      double precision dx,dy,dt,cfl, xlower, ylower, t
      double precision dtdx, dtdy

      double precision work(mwork)

      double precision qold(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      double precision fp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision fm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision gm(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)


c     # Local variables
      integer i0faddm, i0faddp, i0gaddm, i0gaddp
      integer i0q1d, i0dtdx1, i0dtdy1
      integer i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl
      integer i0auxr

      integer i,j,m

c     Needed by Riemann solvers.  This should be fixed later by a 'context'
c     for a Riemann solver.
      double precision dtcom, dxcom,dycom,tcom
      integer icom, jcom
      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      ierror = 0

c     This should be set to actual time, in case the user wants it
c     it for some reason in the Riemann solver.
      tcom = t
      dxcom = dx
      dycom = dy

      maxmx = mx
      maxmy = my

c     # Set up work arrays

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
         ierror = 1
         return
      endif


c     # Include four additional arguments to avoid need for
c     # global array
      if (use_fwaves .eq. 0) then
          call clawpack46_step2(maxm,maxmx,maxmy,meqn,maux, mbc,
     &          mx,my, qold,aux,dx,dy,dt,
     &          cfl,fm,fp,gm,gp,
     &          work(i0faddm),work(i0faddp),
     &          work(i0gaddm),work(i0gaddp),
     &          work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &          work(i0aux1),work(i0aux2),work(i0aux3),
     &          work(i0next),mwork1,rpn2,rpt2,flux2,
     &          mwaves,mcapa,method,mthlim,block_corner_count,ierror)
      else
c         // Only difference is that we pass in rpn2fw, rpt2fw
c         // instead of rpn2 and rpt2 (they have difference signatures). 
c         // But fortunately, Fortran doesn't require separate 
c         // signatures 
          call clawpack46_step2(maxm,maxmx,maxmy,meqn,maux, mbc,
     &          mx,my, qold,aux,dx,dy,dt,
     &          cfl,fm,fp,gm,gp,
     &          work(i0faddm),work(i0faddp),
     &          work(i0gaddm),work(i0gaddp),
     &          work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &          work(i0aux1),work(i0aux2),work(i0aux3),
     &          work(i0next),mwork1,rpn2fw,rpt2fw,flux2,
     &          mwaves,mcapa,method,mthlim,block_corner_count,ierror)
      endif

c     # update q
      dtdx = dt/dx
      dtdy = dt/dy
      do m = 1,meqn
         do i = 1-mbc,mx+mbc-1
            do j = 1-mbc,my+mbc-1
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
