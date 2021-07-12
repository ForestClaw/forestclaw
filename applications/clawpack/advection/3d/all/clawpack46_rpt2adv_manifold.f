      subroutine clawpack46_rpt2adv_manifold(ixy,maxm,meqn,
     &      mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
     &      imp,asdq,bmasdq,bpasdq)
      implicit none

      integer ixy,maxm, meqn, mwaves, mbc, mx, imp, m
      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)

      integer iface, i, i1

      iface = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
         do m = 1,meqn
            bmasdq(i,m) = min(aux2(i1,iface+1), 0.d0) * asdq(i,m)
            bpasdq(i,m) = max(aux3(i1,iface+1), 0.d0) * asdq(i,m)
         enddo
      enddo

      return
      end
