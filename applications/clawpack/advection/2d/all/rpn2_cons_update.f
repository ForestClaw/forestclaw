      subroutine  rpn2_cons_update(meqn,maux, idir,
     &                ql,qr,auxl,auxr,fluxdiff)

      implicit none

      integer meqn,maux,idir, ixy
      double precision ql(meqn), qr(meqn), auxl(maux), auxr(maux)
      double precision fluxdiff(meqn)
      double precision wave, s

c     # Flux function : ur*qr - ul*ql = apdq + amdq

      
      ixy = idir + 1

      wave = qr(1) - ql(1)
      s = auxl(ixy + 1)

      fluxdiff(1) = s*wave

      end
