      subroutine rpn2cons_qs(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
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
      double precision qrr,qll,ur,ul, uhat, umax


      iface = ixy
      do i = 2-mbc, mx+mbc

          ur = auxl(i,iface)
          ul = auxr(i-1,iface)

          qrr = ql(i,1)
          qll = qr(i-1,1)

          umax = max(ur,ul)
          if (umax .ge. 0) then
c             # At least one non-negative velocity            
              if (min(ur,ul) .ge. 0) then
c                 # Both ur, ul >= 0                
                  uhat = ur
              else
c                 # Choose positive speed for convergent/divergent cases                
                  uhat = umax
              endif
              uhat = max(uhat,1d-10)   !! Avoid division by zero
              amdq(i,1) = 0.d0
              apdq(i,1) = ur*qrr - ul*qll
          else
c             # Both ur,ul < 0            
              uhat = ul
              amdq(i,1) = ur*qrr - ul*qll
              apdq(i,1) = 0.d0
          endif
          s(i,1) = uhat  
          wave(i,1,1) = (ur*qrr - ul*qll)/uhat          
      enddo

      return
      end
