      subroutine rpn2cons_fw_manifold(ixy,maxm,meqn,mwaves,mbc,
     &                            mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,maux, mx
      integer ixy

      double precision fwave(1-mbc:maxm+mbc,meqn,mwaves)   
      double precision    s(1-mbc:maxm+mbc,mwaves)
      double precision   ql(1-mbc:maxm+mbc,meqn)
      double precision   qr(1-mbc:maxm+mbc,meqn)
      double precision amdq(1-mbc:maxm+mbc,meqn)
      double precision apdq(1-mbc:maxm+mbc,meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)

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
         g = auxl(i,12 + 2*idir)  

c        # Cell-centered values
         ur(1) = auxl(i,2)
         ur(2) = auxl(i,3)

         ul(1) = auxr(i-1,2)
         ul(2) = auxr(i-1,3)

c        # left edge   : 4,5
c        # bottom-edge : 6,7
         nv(1) = auxl(i,4 + 4*idir)
         nv(2) = auxl(i,5 + 4*idir)

         urrot = g*(nv(1)*ur(1) + nv(2)*ur(2))
         ulrot = g*(nv(1)*ul(1) + nv(2)*ul(2))

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
         fwave(i,1,1) = urrot*qrr - ulrot*qll
         s(i,1) = uhat
      enddo


      return
      end
