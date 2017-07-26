      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &		       wave,s,amdq,apdq)

      implicit double precision (a-h,o-z)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
      dimension auxl(1-mbc:maxmx+mbc, *)
      dimension auxr(1-mbc:maxmx+mbc, *)

      do i = 2-mbc,mx+mbc
         ur = auxl(i,1)
         ul = auxr(i-1,1)

         qrr = ql(i,1)
         qll = qr(i-1,1)

         uhat = ur    !! left edge value

         wave(i,1,1) = qrr-qll
         amdq(i,1) = min(uhat,0.d0)*wave(i,1,1)
         apdq(i,1) = max(uhat,0.d0)*wave(i,1,1)
         s(i,1) = uhat
      enddo
c
c         
c
      return
      end



