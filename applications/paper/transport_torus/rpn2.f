      subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,
     &                mx,ql,qr,
     &                auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,maux, mx
      integer ixy

      double precision wave(meqn,mwaves, 1-mbc:maxm+mbc)   
      double precision    s(mwaves,1-mbc:maxm+mbc)
      double precision   ql(meqn,1-mbc:maxm+mbc)
      double precision   qr(meqn,1-mbc:maxm+mbc)
      double precision amdq(meqn,1-mbc:maxm+mbc)
      double precision apdq(meqn,1-mbc:maxm+mbc)
      double precision auxl(maux,1-mbc:maxm+mbc)
      double precision auxr(maux,1-mbc:maxm+mbc)

      call rpn2cons_fw_manifold(ixy,maxm,meqn,mwaves,maux,mbc,
     &                         mx,ql,qr, auxl,auxr,wave,s,amdq,apdq)

      end