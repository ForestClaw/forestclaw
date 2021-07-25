      subroutine clawpack5_rpn2cons_manifold(ixy,maxm,meqn,
     &         mwaves,maux,mbc, mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer ixy, maxm, meqn, mwaves,maux, mbc, mx

      double precision wave(meqn, mwaves, 1-mbc:maxm+mbc)   
      double precision    s(mwaves, 1-mbc:maxm+mbc)
      double precision   ql(meqn, 1-mbc:maxm+mbc)
      double precision   qr(meqn, 1-mbc:maxm+mbc)
      double precision amdq(meqn, 1-mbc:maxm+mbc)
      double precision apdq(meqn, 1-mbc:maxm+mbc)
      double precision auxl(maux, 1-mbc:maxm+mbc)
      double precision auxr(maux, 1-mbc:maxm+mbc)

      integer :: icom,jcom
      double precision :: dtcom,dxcom,dycom,tcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      integer i, idir
      double precision qll,qrr
      double precision urrot, ulrot, g, uhat

      idir = ixy-1
      do i = 2-mbc, mx+mbc
         g = auxl(6+idir,i)  !! Edge length

c        # left-right : 2,3
c        # top-bottom : 4,5         
         urrot = g*auxl(2 + 2*idir,i)   !! Left edge of right cell
         ulrot = g*auxl(3 + 2*idir,i-1)   !! Right edge of left cell

         qrr = ql(1,i)
         qll = qr(1,i-1)

c        # Use Roe-average values         
         uhat = (ulrot + urrot)/2.d0

         if (uhat .ge. 0) then
            amdq(1,i) = 0.d0
            apdq(1,i) = urrot*qrr - ulrot*qll
         else
            amdq(1,i) = urrot*qrr - ulrot*qll
            apdq(1,i) = 0.d0
         endif
         wave(1,1,i) = urrot*qrr - ulrot*qll
         s(1,i) = uhat
      enddo


      return
      end
