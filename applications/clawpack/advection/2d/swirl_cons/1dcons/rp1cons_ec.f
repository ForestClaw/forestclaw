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
      double precision qrr,qll,ur,ul,uavg,uhat,qhat,fstar


      do i = 2-mbc,mx+mbc
          ur = auxl(i,1)    !! cell-centered.
          ul = auxr(i-1,1)
          uhat = (ur + ul)/2.d0

          qrr = ql(i,1)
          qll = qr(i-1,1)

          wave(i,1,1) = qrr - qll
          s(i,1) = uhat

          if (uhat .ge. 0.d0) then
              fstar = uhat*qll
          else
              fstar = uhat*qrr
          endif

          amdq(i,1) = fstar - ul*qll              !flux
          apdq(i,1) = ur*qrr - fstar            !-flux
      enddo
c
      return
      end



