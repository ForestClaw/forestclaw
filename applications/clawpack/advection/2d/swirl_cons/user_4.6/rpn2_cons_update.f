      subroutine  rpn2_cons_update(meqn,maux, idir,
     &                ql,qr,auxl,auxr,fluxdiff)

      implicit none

      integer meqn,maux,idir, ixy
      double precision ql(meqn), qr(meqn), auxl(maux), auxr(maux)
      double precision fluxdiff(meqn)

c     # Flux function : ur*qr - ul*ql = apdq + amdq

      
      ixy = idir + 1
      fluxdiff(1) = auxr(ixy)*qr(1) - auxl(ixy)*ql(1)


      end
