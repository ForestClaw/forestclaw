      subroutine rpn2cons_cc(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
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


      integer i, iface
      double precision qstar, qrr,qll,ur,ul, uhat, uavg


      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom

      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

C     # We need to impose flux continuity at the "zero-wave"
C     # ur > 0 : us*qs - ul*qll = 0 -->  us = ur; qs = ul*qll/ur;
C     # ur > 0 : us*qs - ul*qll = 0 -->  us = ur; qs = ul*qll/ur;


      iface = ixy
      do i = 2-mbc, mx+mbc

          ur = auxl(i,iface)
          ul = auxr(i-1,iface)

          qrr = ql(i,1)
          qll = qr(i-1,1)

c         # ur and ul have the same sign
          if (ur .gt. 0 .and. ul .gt. 0) then
              wave(i,1,1) = (ur*qrr - ul*qll)/ur
              s(i,1) = ur
              amdq(i,1) = 0.d0
              apdq(i,1) = ur*qrr - ul*qll
          elseif (ur .lt. 0 .and. ul .lt. 0) then
              wave(i,1,1) = (ur*qrr - ul*qll)/ul
              s(i,1) = ul
              amdq(i,1) = ur*qrr - ul*qll
              apdq(i,1) = 0.d0
          else
C             # ul and ur have different signs
              uavg = (ur+ul)/2.d0
              uhat = max(abs(uavg),1d-12)
              s(i,1) = uhat
              wave(i,1,1) = (ur*qrr - ul*qll)/uhat
              amdq(i,1) = -ul*qll
              apdq(i,1) = ur*qrr
          endif

      enddo

      return
      end
