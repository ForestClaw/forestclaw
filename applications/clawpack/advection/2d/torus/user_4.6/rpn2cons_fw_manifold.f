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


      integer i, k, iface, m, idir
      double precision uhat,qll,qrr, uvecl(3),uvecr(3)
      double precision urrot, ulrot, g, nv(3)

      idir = ixy-1     !! 0 for the x-face; 1 for the y-face
      do i = 2-mbc, mx+mbc

c        # Get cell-centered velocities         
         uvecr(1) = auxl(i,4)
         uvecr(2) = auxl(i,5)
         uvecr(3) = auxl(i,6)

         uvecl(1) = auxr(i-1,4)
         uvecl(2) = auxr(i-1,5)
         uvecl(3) = auxr(i-1,6)

c        # x-edge lengths (7)
c        # y-edge lengths (8)
         g = auxl(i,7+idir)

c        # Get scaled edge normals
c        #    --- x-face : (9,10,11)     
c        #    --- y-face : (12,13,14)     
         do k = 1,3
            if (idir .eq. 0) then
               nv(k) = g*auxl(i,9+k-1)
            else
               nv(k) = g*auxl(i,12+k-1)
            endif
         enddo

         urrot = 0
         ulrot = 0
         do k = 1,3
            urrot = urrot + nv(k)*uvecr(k)
            ulrot = ulrot + nv(k)*uvecl(k)
         enddo

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
