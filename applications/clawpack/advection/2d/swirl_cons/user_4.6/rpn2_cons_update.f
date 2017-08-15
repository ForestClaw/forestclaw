      subroutine  rpn2_cons_update(meqn,maux, idir,
     &                ql,qr,auxl,auxr,fluxdiff)

      implicit none

      integer meqn,maux,idir, iface
      double precision ql(meqn), qr(meqn), auxl(maux), auxr(maux)
      double precision fluxdiff(meqn)

c     # Flux function : ur*qr - ul*ql = apdq + amdq

      
      iface = idir + 1
      fluxdiff(1) = auxr(iface)*qr(1) - auxl(iface)*ql(1)


      end
