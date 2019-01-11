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
      double precision uhat,qll,qrr
      double precision urrot, ulrot

      idir = ixy-1     !! 0 for the x-face; 1 for the y-face
      do i = 2-mbc, mx+mbc

c        # Get cell-centered velocities         
c         uvecr(1) = auxl(i,4)
c         uvecr(2) = auxl(i,5)
c         uvecr(3) = auxl(i,6)

c         uvecl(1) = auxr(i-1,4)
c         uvecl(2) = auxr(i-1,5)
c         uvecl(3) = auxr(i-1,6)

c        # x-edge lengths (7)
c        # y-edge lengths (8)
c         g = auxl(i,7+idir)

c        # Get scaled edge normals
c        #    --- x-face : (9,10,11)     
c        #    --- y-face : (12,13,14)     

c         do k = 1,3
c            if (idir .eq. 0) then
c               nv(k) = g*auxl(i,9+k-1)
c            else
c               nv(k) = g*auxl(i,12+k-1)
c            endif
c         enddo

c         urrot = 0
c         ulrot = 0
c         do k = 1,3
c            urrot = urrot + nv(k)*uvecr(k)
c            ulrot = ulrot + nv(k)*uvecl(k)
c         enddo

         urrot = auxl(i,  2 + 2*idir)   !! Left edge of right cell
         ulrot = auxl(i-1,3 + 2*idir)   !! Right edge of left cell

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
