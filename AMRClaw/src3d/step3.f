c
c     ==================================================================
      subroutine step3(maxm,maxmx,maxmy,maxmz,meqn,maux,mbc,mx,my,mz,
     &		       qold,aux,dx,dy,dz,dt,cflgrid,
     &                 fm,fp,gm,gp,hm,hp,
     &                 faddm,faddp,gadd,hadd,
     &                 q1d,dtdx1d,dtdy1d,dtdz1d,
     &                 aux1,aux2,aux3,work,mwork,rpn3,rpt3,rptt3)
c     ==================================================================

c     # clawpack routine ...  modified for AMRCLAW

c
c     # Take one time step, updating q.
c     # On entry, qold gives
c     #    initial data for this step
c     #    and is unchanged in this version.
c
c     # fm, fp are fluxes to left and right of single cell edge
c
c     # See the flux3 documentation for more information.
c
c
      implicit none
      external rpn3,rpt3,rptt3
c      include  "call.i"
      include "claw.i"


      integer maxm, maxmx, maxmy, maxmz, meqn, maux, mbc, mx,
     &      my, mz, mwork

      double precision dx, dy, dz, dt, cflgrid

      double precision qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      double precision  q1d(1-mbc:maxm+mbc, meqn)
c
      double precision   fm(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision   fp(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision   gm(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision   gp(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision   hm(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision   hp(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, meqn)
c
      double precision faddm(1-mbc:maxm+mbc, meqn)
      double precision faddp(1-mbc:maxm+mbc, meqn)
      double precision  gadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
      double precision  hadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
c
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &              1-mbc:maxmz+mbc, *)
      double precision aux1(1-mbc:maxm+mbc, maux, 3)
      double precision aux2(1-mbc:maxm+mbc, maux, 3)
      double precision aux3(1-mbc:maxm+mbc, maux, 3)
      double precision dtdx1d(1-mbc:maxm+mbc)
      double precision dtdy1d(1-mbc:maxm+mbc)
      double precision dtdz1d(1-mbc:maxm+mbc)
      double precision work(mwork)


      integer i0wave,  i0s, i0amdq, i0apdq, i0cqxx, i0bmamdq,  i0bmapdq,
     &      i0bpamdq, i0bpapdq, i0cmamdq, i0cmapdq,i0cpamdq, i0cpapdq,
     &      i0cmamdq2, i0cmapdq2, i0cpamdq2,i0cpapdq2, i0bmcqxxm,
     &      i0bmcqxxp, i0bpcqxxm, i0bpcqxxp,i0cmcqxxm, i0cmcqxxp,
     &      i0cpcqxxm,  i0cpcqxxp, i0bmcmamdq,i0bmcmapdq, i0bpcmamdq,
     &      i0bpcmapdq, i0bmcpamdq, i0bmcpapdq,i0bpcpamdq,i0bpcpapdq,
     &      iused

      double precision dtdx, dtdy, dtdz, cfl1d
      integer m,i,j,k, ma

      double precision dtcom, dxcom, dycom, dzcom, tcom
      integer icom, jcom, kcom



      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
c
c     # store mesh parameters that may be needed in Riemann solver but not
c     # passed in...
      dxcom = dx
      dycom = dy
      dtcom = dt

c
c
c     # partition work array into pieces needed for local storage in
c     # flux3 routine.  Find starting index of each piece:
c
      i0wave     = 1
      i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
      i0amdq     = i0s        + (maxm+2*mbc)*mwaves
      i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
      i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
      i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
      i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
      i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
      i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
      i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
      i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
      i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
      i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn

      i0bmcqxxm   = i0cpapdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxp   = i0bmcqxxm  + (maxm+2*mbc)*meqn
      i0bpcqxxm   = i0bmcqxxp  + (maxm+2*mbc)*meqn
      i0bpcqxxp   = i0bpcqxxm  + (maxm+2*mbc)*meqn
      i0cmcqxxm   = i0bpcqxxp  + (maxm+2*mbc)*meqn
      i0cmcqxxp   = i0cmcqxxm  + (maxm+2*mbc)*meqn
      i0cpcqxxm   = i0cmcqxxp  + (maxm+2*mbc)*meqn
      i0cpcqxxp   = i0cpcqxxm  + (maxm+2*mbc)*meqn

      i0bmcmamdq = i0cpcqxxp  + (maxm+2*mbc)*meqn
      i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
      i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
      i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
      i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
      i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
      i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
      i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
      iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen if mwork is set properly in stepgrid3
	 write(6,*) '*** not enough work space in step3'
	 write(6,*) '*** check parameter mwork in stepgrid3'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
	 stop
      endif
c
c
      cflgrid = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
c
      do 10 m=1,meqn
         do 10 i=1-mbc,mx+mbc
            do 10 j=1-mbc,my+mbc
              do 10 k=1-mbc,mz+mbc
               fm(i,j,k,m) = 0.d0
               fp(i,j,k,m) = 0.d0
               gm(i,j,k,m) = 0.d0
               gp(i,j,k,m) = 0.d0
               hm(i,j,k,m) = 0.d0
               hp(i,j,k,m) = 0.d0
   10          continue

c
      if (mcapa.eq.0) then
c        # no capa array:
	 do 5 i=1-mbc,maxm+mbc
	    dtdx1d(i) = dtdx
	    dtdy1d(i) = dtdy
	    dtdz1d(i) = dtdz
    5       continue
	 endif
c
c
c     # perform x-sweeps
c     ==================
c
      do 50 k = 0,mz+1
         do 50 j = 0,my+1
c
 	    do 20 m=1,meqn
 	       do 20 i = 1-mbc, mx+mbc
c                 # copy data along a slice into 1d array:
 	          q1d(i,m) = qold(i,j,k,m)
   20          continue
c
         if (mcapa.gt.0)  then
           do 23 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,k,mcapa)
   23      continue
         endif
c
         if (maux .gt. 0)  then
             do 22 ma=1,maux
               do 22 i = 1-mbc, mx+mbc
                 aux1(i,ma,1) = aux(i,j-1,k-1,ma)
                 aux1(i,ma,2) = aux(i,j-1,k,ma)
                 aux1(i,ma,3) = aux(i,j-1,k+1,ma)
                 aux2(i,ma,1) = aux(i,j,k-1,ma)
                 aux2(i,ma,2) = aux(i,j,k,ma)
                 aux2(i,ma,3) = aux(i,j,k+1,ma)
                 aux3(i,ma,1) = aux(i,j+1,k-1,ma)
                 aux3(i,ma,2) = aux(i,j+1,k,ma)
                 aux3(i,ma,3) = aux(i,j+1,k+1,ma)
   22          continue
           endif
c
c	    # Store the value of j and k along this slice in the common block
c           # comxyt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
	    jcom = j
            kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along
c           # this slice:
c
            call flux3(1,maxm,meqn,maux,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
	    cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
         do 25 m=1,meqn
            do 25 i=1,mx+1
               fm(i,j,k,m) = fm(i,j,k,m) + faddm(i,m)
               fp(i,j,k,m) = fp(i,j,k,m) + faddp(i,m)
c
               gm(i,j  ,k-1,m) = gm(i,j  ,k-1,m) + gadd(i,m,1,-1)
               gp(i,j  ,k-1,m) = gp(i,j  ,k-1,m) + gadd(i,m,1,-1)
               gm(i,j  ,k,  m) = gm(i,j  ,k,  m) + gadd(i,m,1, 0)
               gp(i,j  ,k,  m) = gp(i,j  ,k,  m) + gadd(i,m,1, 0)
               gm(i,j  ,k+1,m) = gm(i,j  ,k+1,m) + gadd(i,m,1, 1)
               gp(i,j  ,k+1,m) = gp(i,j  ,k+1,m) + gadd(i,m,1, 1)
c
               gm(i,j+1,k-1,m) = gm(i,j+1,k-1,m) + gadd(i,m,2,-1)
               gp(i,j+1,k-1,m) = gp(i,j+1,k-1,m) + gadd(i,m,2,-1)
               gm(i,j+1,k,  m) = gm(i,j+1,k,  m) + gadd(i,m,2, 0)
               gp(i,j+1,k,  m) = gp(i,j+1,k,  m) + gadd(i,m,2, 0)
               gm(i,j+1,k+1,m) = gm(i,j+1,k+1,m) + gadd(i,m,2, 1)
               gp(i,j+1,k+1,m) = gp(i,j+1,k+1,m) + gadd(i,m,2, 1)
c
               hm(i,j-1,k  ,m) = hm(i,j-1,k  ,m) + hadd(i,m,1,-1)
               hp(i,j-1,k  ,m) = hp(i,j-1,k  ,m) + hadd(i,m,1,-1)
               hm(i,j  ,k  ,m) = hm(i,j  ,k  ,m) + hadd(i,m,1, 0)
               hp(i,j  ,k  ,m) = hp(i,j  ,k  ,m) + hadd(i,m,1, 0)
               hm(i,j+1,k  ,m) = hm(i,j+1,k  ,m) + hadd(i,m,1, 1)
               hp(i,j+1,k  ,m) = hp(i,j+1,k  ,m) + hadd(i,m,1, 1)
c
               hm(i,j-1,k+1,m) = hm(i,j-1,k+1,m) + hadd(i,m,2,-1)
               hp(i,j-1,k+1,m) = hp(i,j-1,k+1,m) + hadd(i,m,2,-1)
               hm(i,j  ,k+1,m) = hm(i,j  ,k+1,m) + hadd(i,m,2, 0)
               hp(i,j  ,k+1,m) = hp(i,j  ,k+1,m) + hadd(i,m,2, 0)
               hm(i,j+1,k+1,m) = hm(i,j+1,k+1,m) + hadd(i,m,2, 1)
               hp(i,j+1,k+1,m) = hp(i,j+1,k+1,m) + hadd(i,m,2, 1)

   25          continue
   50    continue
c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 k = 0, mz+1
         do 100 i = 0, mx+1
c
	    do 70 m=1,meqn
	       do 70 j = 1-mbc, my+mbc
c                 # copy data along a slice into 1d array:
	          q1d(j,m) = qold(i,j,k,m)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,k,mcapa)
   71      continue
         endif
c
         if (maux .gt. 0)  then
             do 72 ma=1,maux
               do 72 j = 1-mbc, my+mbc
                 aux1(j,ma,1) = aux(i-1,j,k-1,ma)
                 aux1(j,ma,2) = aux(i,j,k-1,ma)
                 aux1(j,ma,3) = aux(i+1,j,k-1,ma)
                 aux2(j,ma,1) = aux(i-1,j,k,ma)
                 aux2(j,ma,2) = aux(i,j,k,ma)
                 aux2(j,ma,3) = aux(i+1,j,k,ma)
                 aux3(j,ma,1) = aux(i-1,j,k+1,ma)
                 aux3(j,ma,2) = aux(i,j,k+1,ma)
                 aux3(j,ma,3) = aux(i+1,j,k+1,ma)
   72          continue
         endif
c
c	    # Store the value of i and k along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(2,maxm,meqn,maux,mbc,my,
     &                 q1d,dtdy1d,dtdz,dtdx,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the g-fluxes
c           # gadd - modifies the h-fluxes
c           # hadd - modifies the f-fluxes
c
         do 75 m=1,meqn
            do 75 j=1,my+1
               gm(i,j,k,m) = gm(i,j,k,m) + faddm(j,m)
               gp(i,j,k,m) = gp(i,j,k,m) + faddp(j,m)
c
               hm(i-1,j,k  ,m) = hm(i-1,j,k  ,m) + gadd(j,m,1,-1)
               hp(i-1,j,k  ,m) = hp(i-1,j,k  ,m) + gadd(j,m,1,-1)
               hm(i  ,j,k  ,m) = hm(i  ,j,k  ,m) + gadd(j,m,1, 0)
               hp(i  ,j,k  ,m) = hp(i  ,j,k  ,m) + gadd(j,m,1, 0)
               hm(i+1,j,k  ,m) = hm(i+1,j,k  ,m) + gadd(j,m,1, 1)
               hp(i+1,j,k  ,m) = hp(i+1,j,k  ,m) + gadd(j,m,1, 1)
c
               hm(i-1,j,k+1,m) = hm(i-1,j,k+1,m) + gadd(j,m,2,-1)
               hp(i-1,j,k+1,m) = hp(i-1,j,k+1,m) + gadd(j,m,2,-1)
               hm(i  ,j,k+1,m) = hm(i  ,j,k+1,m) + gadd(j,m,2, 0)
               hp(i  ,j,k+1,m) = hp(i  ,j,k+1,m) + gadd(j,m,2, 0)
               hm(i+1,j,k+1,m) = hm(i+1,j,k+1,m) + gadd(j,m,2, 1)
               hp(i+1,j,k+1,m) = hp(i+1,j,k+1,m) + gadd(j,m,2, 1)
c
               fm(i  ,j,k-1,m) = fm(i  ,j,k-1,m) + hadd(j,m,1,-1)
               fp(i  ,j,k-1,m) = fp(i  ,j,k-1,m) + hadd(j,m,1,-1)
               fm(i  ,j,k  ,m) = fm(i  ,j,k  ,m) + hadd(j,m,1, 0)
               fp(i  ,j,k  ,m) = fp(i  ,j,k  ,m) + hadd(j,m,1, 0)
               fm(i  ,j,k+1,m) = fm(i  ,j,k+1,m) + hadd(j,m,1, 1)
               fp(i  ,j,k+1,m) = fp(i  ,j,k+1,m) + hadd(j,m,1, 1)
c
               fm(i+1,j,k-1,m) = fm(i+1,j,k-1,m) + hadd(j,m,2,-1)
               fp(i+1,j,k-1,m) = fp(i+1,j,k-1,m) + hadd(j,m,2,-1)
               fm(i+1,j,k  ,m) = fm(i+1,j,k  ,m) + hadd(j,m,2, 0)
               fp(i+1,j,k  ,m) = fp(i+1,j,k  ,m) + hadd(j,m,2, 0)
               fm(i+1,j,k+1,m) = fm(i+1,j,k+1,m) + hadd(j,m,2, 1)
               fp(i+1,j,k+1,m) = fp(i+1,j,k+1,m) + hadd(j,m,2, 1)
c

   75          continue
c
  100    continue
c
c
c
c     # perform z sweeps
c     ==================
c
c
      do 150 j = 0, my+1
         do 150 i = 0, mx+1
c
	    do 110 m=1,meqn
	       do 110 k = 1-mbc, mz+mbc
c                 # copy data along a slice into 1d array:
	          q1d(k,m) = qold(i,j,k,m)
 110           continue
c
         if (mcapa.gt.0)  then
           do 130 k = 1-mbc, mz+mbc
               dtdz1d(k) = dtdz / aux(i,j,k,mcapa)
 130       continue
         endif
c
         if (maux .gt. 0)  then
             do 131 ma=1,maux
               do 131 k = 1-mbc, mz+mbc
                 aux1(k,ma,1) = aux(i-1,j-1,k,ma)
                 aux1(k,ma,2) = aux(i-1,j,k,ma)
                 aux1(k,ma,3) = aux(i-1,j+1,k,ma)
                 aux2(k,ma,1) = aux(i,j-1,k,ma)
                 aux2(k,ma,2) = aux(i,j,k,ma)
                 aux2(k,ma,3) = aux(i,j+1,k,ma)
                 aux3(k,ma,1) = aux(i+1,j-1,k,ma)
                 aux3(k,ma,2) = aux(i+1,j,k,ma)
                 aux3(k,ma,3) = aux(i+1,j+1,k,ma)
  131          continue
           endif
c
c	    # Store the value of i and j along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            jcom = j
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(3,maxm,meqn,maux,mbc,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
	    cflgrid = dmax1(cflgrid,cfl1d)
c
c
c        # update fluxes for use in AMR:
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the h-fluxes
c           # gadd - modifies the f-fluxes
c           # hadd - modifies the g-fluxes
c
         do 125 m=1,meqn
            do 125 k=1,mz+1
               hm(i,j,k,m) = hm(i,j,k,m) + faddm(k,m)
               hp(i,j,k,m) = hp(i,j,k,m) + faddp(k,m)
c
               fm(i  ,j-1,k,m) = fm(i  ,j-1,k,m) + gadd(k,m,1,-1)
               fp(i  ,j-1,k,m) = fp(i  ,j-1,k,m) + gadd(k,m,1,-1)
               fm(i  ,j  ,k,m) = fm(i  ,j  ,k,m) + gadd(k,m,1, 0)
               fp(i  ,j  ,k,m) = fp(i  ,j  ,k,m) + gadd(k,m,1, 0)
               fm(i  ,j+1,k,m) = fm(i  ,j+1,k,m) + gadd(k,m,1, 1)
               fp(i  ,j+1,k,m) = fp(i  ,j+1,k,m) + gadd(k,m,1, 1)
c
               fm(i+1,j-1,k,m) = fm(i+1,j-1,k,m) + gadd(k,m,2,-1)
               fp(i+1,j-1,k,m) = fp(i+1,j-1,k,m) + gadd(k,m,2,-1)
               fm(i+1,j  ,k,m) = fm(i+1,j  ,k,m) + gadd(k,m,2, 0)
               fp(i+1,j  ,k,m) = fp(i+1,j  ,k,m) + gadd(k,m,2, 0)
               fm(i+1,j+1,k,m) = fm(i+1,j+1,k,m) + gadd(k,m,2, 1)
               fp(i+1,j+1,k,m) = fp(i+1,j+1,k,m) + gadd(k,m,2, 1)
c
               gm(i-1,j  ,k,m) = gm(i-1,j  ,k,m) + hadd(k,m,1,-1)
               gp(i-1,j  ,k,m) = gp(i-1,j  ,k,m) + hadd(k,m,1,-1)
               gm(i  ,j  ,k,m) = gm(i  ,j  ,k,m) + hadd(k,m,1, 0)
               gp(i  ,j  ,k,m) = gp(i  ,j  ,k,m) + hadd(k,m,1, 0)
               gm(i+1,j  ,k,m) = gm(i+1,j  ,k,m) + hadd(k,m,1, 1)
               gp(i+1,j  ,k,m) = gp(i+1,j  ,k,m) + hadd(k,m,1, 1)
c
               gm(i-1,j+1,k,m) = gm(i-1,j+1,k,m) + hadd(k,m,2,-1)
               gp(i-1,j+1,k,m) = gp(i-1,j+1,k,m) + hadd(k,m,2,-1)
               gm(i  ,j+1,k,m) = gm(i  ,j+1,k,m) + hadd(k,m,2, 0)
               gp(i  ,j+1,k,m) = gp(i  ,j+1,k,m) + hadd(k,m,2, 0)
               gm(i+1,j+1,k,m) = gm(i+1,j+1,k,m) + hadd(k,m,2, 1)
               gp(i+1,j+1,k,m) = gp(i+1,j+1,k,m) + hadd(k,m,2, 1)
c
  125          continue
  150    continue
c
c
      return
      end
