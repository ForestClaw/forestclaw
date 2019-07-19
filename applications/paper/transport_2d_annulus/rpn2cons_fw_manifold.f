      subroutine rpn2cons_fw_manifold(ixy,maxm,meqn,mwaves,maux,mbc,
     &                            mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,maux, mx
      integer ixy

      double precision fwave(meqn,mwaves, 1-mbc:maxm+mbc)   
      double precision    s(mwaves,1-mbc:maxm+mbc)
      double precision   ql(meqn,1-mbc:maxm+mbc)
      double precision   qr(meqn,1-mbc:maxm+mbc)
      double precision amdq(meqn,1-mbc:maxm+mbc)
      double precision apdq(meqn,1-mbc:maxm+mbc)
      double precision auxl(maux,1-mbc:maxm+mbc)
      double precision auxr(maux,1-mbc:maxm+mbc)

      integer i, iface, m, idir
      double precision qll,qrr
      double precision urrot, ulrot, g, uhat
      double precision nv(2), ur(2), ul(2)
      double precision g1, g2, n1, n2

      logical qad_debug
      common /debug_common/ qad_debug

      idir = ixy-1
      do i = 2-mbc, mx+2
         !! Edge length;  assumes that edge length is stored at the 
         !! left edge.
         g = auxl(12 + 2*idir,i)  

c        # Cell-centered values
         ur(1) = auxl(2,i)
         ur(2) = auxl(3,i)

         ul(1) = auxr(2,i-1)
         ul(2) = auxr(3,i-1)

c        # left edge   : 4,5
c        # bottom-edge : 6,7
         nv(1) = auxl(4 + 4*idir,i)
         nv(2) = auxl(5 + 4*idir,i)

         urrot = g*(nv(1)*ur(1) + nv(2)*ur(2))
         ulrot = g*(nv(1)*ul(1) + nv(2)*ul(2))

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
         fwave(1,1,i) = urrot*qrr - ulrot*qll
         s(1,i) = uhat
      enddo


      return
      end
