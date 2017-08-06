      subroutine  rpn2_cons_update(meqn,maux, idir,
     &                ql,qr,auxl,auxr,fluxdiff)

      implicit none

      integer meqn,maux,idir
      double precision ql(meqn), qr(meqn), auxl(maux), auxr(maux)
      double precision fluxdiff(meqn)

c     # Flux function : ur*qr - ul*ql = apdq + amdq

      
      if (idir .eq. 0) then
          fluxdiff(1) = auxr(1)*qr(1) - auxl(1)*ql(1)
      else
          fluxdiff(1) = auxr(2)*qr(1) - auxl(2)*ql(1)
      endif


      end
