c      subroutine rpt2adv(ixy,maxm,meqn,
c     &      mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
c     &      imp,asdq,bmasdq,bpasdq)
      subroutine rpt2adv_manifold(ixy,imp,maxm,meqn,mwaves,maux,
     &               mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      implicit none

      integer ixy,maxm, maux, meqn, mwaves, mbc, mx, imp, m
      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)

      integer iface, i, i1
      double precision gp, gm, vm, vp

      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

         !! Note that edge lengths are stored in 6,7,8,9 (different from
         !! ForestClaw).  The reason for this is clean in rpn2_qad.f
         if (ixy .eq. 1) then
            gm = aux2(8,i1)
            gp = aux3(8,i1)
            vm = gm*aux2(4,i1)   !! cell-velocity projected to bottom edge
            vp = gp*aux2(5,i1)   !! cell-velocity projected to top edge
         else
            gm = aux2(6,i1)
            gp = aux3(6,i1)
            vm = gm*aux2(2,i1)   !! cell-velocity projected to left edge
            vp = gp*aux2(3,i1)   !! cell-velocity projected to right edge
         endif
         do m = 1,meqn
            bmasdq(m,i) = min(vm,0.d0)*asdq(m,i)
            bpasdq(m,i) = max(vp,0.d0)*asdq(m,i)
         enddo
      enddo

      return
      end
