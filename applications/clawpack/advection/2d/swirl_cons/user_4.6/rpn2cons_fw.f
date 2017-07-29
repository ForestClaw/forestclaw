      subroutine rpn2cons_fw(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
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


      integer i, iface, m
      double precision uhat,qll,qrr,ul,ur, tol


      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom

      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      iface = ixy
      do i = 2-mbc, mx+mbc

         ur = auxl(i,iface)
         ul = auxr(i-1,iface)

         qrr = ql(i,1)
         qll = qr(i-1,1)

         uhat = (ul + ur)/2.d0

         wave(i,1,1) = ur*qrr - ul*qll
         s(i,1) = uhat
         if (uhat .eq. 0) then
            amdq(i,1) = -ul*qll
            apdq(i,1) = ur*qrr            
         elseif (uhat .lt. 0) then
            amdq(i,1) = wave(i,1,1)
            apdq(i,1) = 0.d0
         else
            amdq(i,1) = 0.d0
            apdq(i,1) = wave(i,1,1)
         endif
      enddo


      return
      end
