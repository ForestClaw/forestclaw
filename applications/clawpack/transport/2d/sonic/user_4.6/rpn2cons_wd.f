      subroutine rpn2cons_wd(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                      auxl,auxr,wave,s,amdq,apdq)

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
      double precision ur,ul,qrr,qll, uhat,uavg,tol

      if ((mwaves .ne. 1) .or. (meqn .ne. 1)) then
         write(6,*) 'rpn2cons : mwaves .ne. 1 or meqn .ne. 1'
         stop
      endif


      iface = ixy
      do i = 2-mbc, mx+mbc
          ur = auxl(i,iface)
          ul = auxr(i-1,iface)

          qrr = ql(i,1)
          qll = qr(i-1,1)

          uavg = (ul + ur)/2.d0

c         # -eps < uavg < eps --> uhat = eps
          tol = 1d-12
          uhat = sign(1.d0,uavg)*max(abs(uavg),tol)

c         # Use of .ge. here is consistent with choice of sign above
          if (uhat .ge. 0) then
              amdq(i,1) = 0
              apdq(i,1) = ur*qrr - ul*qll
          else
              amdq(i,1) = ur*qrr - ul*qll
              apdq(i,1) = 0
          endif
          wave(i,1,1) = (ur*qrr - ul*qll)/uhat   !! Same as above
          s(i,1) = uhat      

c          qhat = (qrr + qll)/2.d0
c
c          uavg = (ur+ul)/2.d0
c          if (uavg .eq. 0) then
c              uhat = 1.d-12
c          else
c              uhat = uavg
c          endif
c
c          s(i,1) = uavg
c          wave(i,1,1) = (ur*qrr - ul*qll)/uhat   !! Same as above
c
c          if (uavg .eq. 0) then
cc             # This definition is needed to guarantee conservation, i.e. 
cc             #       f(qr) - f(ql) = apdq + amdq                      
cc             # and a consistent definition of f(qr), f(ql).  
cc             # If we base apdq,amdq on uhat, we will have inconsistent
cc             # definitions of f(q). 
c              amdq(i,1) = -ul*qhat
c              apdq(i,1) = ur*qhat
c          else
c              amdq(i,1) = min(uhat,0.d0)*wave(i,1,1)
c              apdq(i,1) = max(uhat,0.d0)*wave(i,1,1)
c         endif
      enddo

      return
      end
