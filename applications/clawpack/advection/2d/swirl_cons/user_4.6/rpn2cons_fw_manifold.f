      subroutine rpn2cons_fw_manifold(ixy,maxm,meqn,mwaves,mbc,
     &         mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
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


      integer i, iface, m, idir
      double precision qll,qrr
      double precision urrot, ulrot, g, uhat

      idir = ixy-1
      do i = 2-mbc, mx+mbc
         g = auxl(i,6+idir)  !! Edge length

c        # left-right : 2,3
c        # top-bottom : 4,5         
         urrot = g*auxl(i,  2 + 2*idir)   !! Left edge of right cell
         ulrot = g*auxl(i-1,3 + 2*idir)   !! Right edge of left cell


c         ur = auxl(i,2)
c         ul = auxr(i-1,1)
c         vr = auxl(i,3)
c         vl = auxr(i-1,2)

c        # x-edgelengths (4)
c        # y-edgelengths (5)

c         g = auxl(i,3+iface)

c        # x-face normal : (6,7)
c        # y-face normal : (8,9)      
c         if (ixy .eq. 1) then
c            urrot = g*(auxl(i,6)*ur + auxl(i,7)*vr)
c            ulrot = g*(auxl(i,6)*ul + auxl(i,7)*vl)
c         else
c            urrot = g*(auxl(i,8)*ur + auxl(i,9)*vr)
c            ulrot = g*(auxl(i,8)*ul + auxl(i,9)*vl)
c         endif

         qrr = ql(i,1)
         qll = qr(i-1,1)

c        # Use Roe-average values         
         uhat = (ulrot + urrot)/2.d0

         if (uhat .ge. 0) then
            amdq(i,1) = 0.d0
            apdq(i,1) = urrot*qrr - ulrot*qll
         else
            amdq(i,1) = urrot*qrr - ulrot*qll
            apdq(i,1) = 0.d0
         endif
         wave(i,1,1) = urrot*qrr - ulrot*qll
         s(i,1) = uhat

      enddo


      return
      end
