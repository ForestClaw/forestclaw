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

          uhat = (ul + ur)/2.d0

          wave(i,1,1) = ur*qrr - ul*qll
          s(i,1) = uhat
          if (uhat .eq. 0) then
              amdq(i,1) = -ul*qll
              apdq(i,1) = ur*qrr            
          elseif (uhat .lt. 0) then
              amdq(i,1) = wave(i,1,1)
              apdq(i,1) = 0.d0
          else
              amdq(i,1) = 0.d0
              apdq(i,1) = wave(i,1,1)
          endif

      enddo
c
      return
      end



