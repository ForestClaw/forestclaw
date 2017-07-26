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
      double precision qrr,qll,ur,ul,uavg,uhat


      do i = 2-mbc,mx+mbc
         ur = auxl(i,1)
         ul = auxr(i-1,1)

         qrr = ql(i,1)
         qll = qr(i-1,1)

c         # ur and ul have the same sign
          if (ur .gt. 0 .and. ul .gt. 0) then
              wave(i,1,1) = (ur*qrr - ul*qll)/ur
              s(i,1) = ur
              amdq(i,1) = 0.d0
              apdq(i,1) = ur*qrr - ul*qll
          elseif (ur .lt. 0 .and. ul .lt. 0) then
              wave(i,1,1) = (ur*qrr - ul*qll)/ul
              s(i,1) = ul
              amdq(i,1) = ur*qrr - ul*qll
              apdq(i,1) = 0.d0
          else
C             # ul and ur have different signs
              uavg = (ur+ul)/2.d0
              uhat = max(abs(uavg),1d-12)
              s(i,1) = uhat
              wave(i,1,1) = (ur*qrr - ul*qll)/uhat
              amdq(i,1) = -ul*qll
              apdq(i,1) = ur*qrr
          endif

      enddo
c
c         
c
      return
      end



