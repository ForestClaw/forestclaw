      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &		       wave,s,amdq,apdq)

      implicit none

      integer maxmx, meqn, mwaves, mbc,mx

      double precision   ql(1-mbc:maxmx+mbc, meqn)
      double precision   qr(1-mbc:maxmx+mbc, meqn)
      double precision    s(1-mbc:maxmx+mbc, mwaves)
      double precision wave(1-mbc:maxmx+mbc, meqn, mwaves)
      double precision amdq(1-mbc:maxmx+mbc, meqn)
      double precision apdq(1-mbc:maxmx+mbc, meqn)
      double precision auxl(1-mbc:maxmx+mbc, *)
      double precision auxr(1-mbc:maxmx+mbc, *)

      integer i
      double precision qrr,qll,ur,ul,uavg,uhat,qhat


      do i = 2-mbc,mx+mbc
          ur = auxl(i,1)
          ul = auxr(i-1,1)

          qrr = ql(i,1)
          qll = qr(i-1,1)

          qhat = (qrr + qll)/2.d0

          uavg = (ur+ul)/2.d0
          if (uavg .eq. 0) then
              uhat = 1.d-12
          else
              uhat = uavg
          endif

          s(i,1) = uavg
          wave(i,1,1) = (ur*qrr - ul*qll)/uhat   !! Same as above

          if (uavg .eq. 0) then
c             # This definition is needed to guarantee conservation, i.e. 
c             #       f(qr) - f(ql) = apdq + amdq                      
c             # and a consistent definition of f(qr), f(ql).  
c             # If we base apdq,amdq on uhat, we will have inconsistent
c             # definitions of f(q). 
              amdq(i,1) = -ul*qhat
              apdq(i,1) = ur*qhat
          else
              amdq(i,1) = min(uhat,0.d0)*wave(i,1,1)
              apdq(i,1) = max(uhat,0.d0)*wave(i,1,1)
          endif
      enddo
c
      return
      end



