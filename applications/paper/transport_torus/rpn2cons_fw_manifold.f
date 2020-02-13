      subroutine rpn2cons_fw_manifold(ixy,maxm,meqn,mwaves,maux,mbc,
     &                            mx,ql,qr,
     &                            auxl,auxr,fwave,s,amdq,apdq)
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

c     # Riemann solver assumes that there is unique edge between
c     # each cell.  This is not true when solving between coarse 
c     # and fine grids (in qad fix). 

      idir = ixy-1
      do i = 2-mbc, mx+mbc
         !! Single edge between left and right states
         g = auxl(6 + 2*idir,i)  

c        # Projected velocity stored in each cell.  
c        # 2-3 : velocities project to left-right edges of each cell
c        # 4-5 : velocities project to top-bottom edges of each cell
         urrot = g*auxl(2 + 2*idir,i)     !! Left edge of right cell
         ulrot = g*auxr(3 + 2*idir,i-1)   !! Right edge of left cell

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
