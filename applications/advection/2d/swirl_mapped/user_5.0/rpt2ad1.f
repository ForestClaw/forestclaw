c     =====================================================
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the advection equation.
c

      integer ixy, maux, maxm, meqn, mwaves, mbc, mx, imp, m
      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)

      integer kv, i, i1

      kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
         do m = 1,meqn
            bmasdq(m,i) = min(aux2(kv,i1), 0.d0) * asdq(m,i)
            bpasdq(m,i) = max(aux3(kv,i1), 0.d0) * asdq(m,i)
         enddo
      enddo

      return
      end
