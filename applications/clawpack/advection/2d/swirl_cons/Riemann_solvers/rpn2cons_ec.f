      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &      auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,mx
      integer ixy

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)


      integer i, iface
      double precision f0, vn


      double precision dtcom, dxcom, dycom, dzcom, tcom
      integer icom, jcom, kcom
      integer mq, mw
      common/comxyzt/dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

      integer vflag, get_vflag

      vflag = get_vflag()
      if (vflag .gt. 2) then
         write(6,*) 'rpn2cons_ec : This Riemann solver requires ',
     &         'edge-based velocities.  Set vflag = 1 or 2. '
         write(6,*) 'vflag = ', vflag
         stop
      endif

      iface = ixy
      do i = 2-mbc, mx+mbc
         do mq = 1,meqn
            do mw = 1,meqn
               wave(i,mq,mw) = 0
            enddo
            wave(i,mq,mq) = ql(i,mq) - qr(i-1,mq)
            vn = auxl(i,1 + iface)
            s(i,mq) = vn

            if (vn .ge. 0.d0) then
               f0 = vn*ql(i-1,mq)
            else
               f0 = vn*ql(i,mq)
            endif

            amdq(i,mq) = f0              !flux
            apdq(i,mq) = -f0             !-flux
         enddo
c         vn = auxl(i,10+iface)
c         s(i,1) = vn

c        # f0 is the interface flux
c        # fl, fr are the fluxes at ql, qr
c        # Set f0=0, since it should cancel anyway.

c         if (vn .ge. 0.d0) then
c           f0 = vn*ql(i-1,1)
c         else
c           f0 = vn*ql(i,1)
c         endif
c
c	 amdq(i,1) = f0    !flux
c	 apdq(i,1) = -f0   !-flux
      enddo


      return
      end
