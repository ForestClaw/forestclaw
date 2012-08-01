c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the advection equation.
c

      integer ixy,maxm, meqn, mwaves, mbc, mx, imp
      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, 2)
      double precision   aux2(1-mbc:maxm+mbc, 2)
      double precision   aux3(1-mbc:maxm+mbc, 2)

      integer kv, i, i1

      kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
         bmasdq(i,1) = dmin1(aux2(i1,kv), 0.d0) * asdq(i,1)
         bpasdq(i,1) = dmax1(aux3(i1,kv), 0.d0) * asdq(i,1)
      enddo

      return
      end
