      subroutine rpn2cons(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                       auxl,auxr,wave,s,amdq,apdq)

      implicit none

      integer maxm, mbc,mwaves,meqn,maux, mx
      integer ixy

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)


      integer i, iface, m, mw, version
      double precision ur,ul,qrr,qll, uhat,qhat, tau,uavg
      double precision a,b,c


      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom

      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      if ((mwaves .ne. 1) .or. (meqn .ne. 1)) then
         write(6,*) 'rpn2cons : mwaves .ne. 1 or meqn .ne. 1'
         stop
      endif


      iface = ixy
      do i = 2-mbc, mx+mbc

         ur = auxl(i,iface)
         ul = auxr(i-1,iface)

         qrr = ql(i,1)
         qll = qr(i-1,1)

         qhat = (qrr + qll)/2.d0

         uavg = (ur+ul)/2.d0

         if (uavg .eq. 0) then
            tau = 1e-12
            uhat = sign(1.d0,uavg)*max(abs(uavg),tau) 
         else
            uhat = uavg
         endif
         wave(i,1,1) = qrr-qll + (ur-ul)*qhat/uhat         
         s(i,1) = uhat

         if (uavg .eq. 0) then
            amdq(i,1) = -ul*qhat
            apdq(i,1) = ur*qhat
         else
            amdq(i,1) = min(uhat,0.d0)*wave(i,1,1)
            apdq(i,1) = max(uhat,0.d0)*wave(i,1,1)
         endif
      enddo

      return
      end
