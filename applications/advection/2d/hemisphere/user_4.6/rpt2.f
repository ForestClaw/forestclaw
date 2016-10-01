      subroutine clawpack46_rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &      ql,qr,aux1,aux2,aux3,
     &      imp,asdq,bmasdq,bpasdq)
      implicit none


      integer ixy,imp,maxm,meqn,mwaves,mbc,mx

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)


      double precision sm,sp
      integer i,i1,iface, mq

      iface = 3-ixy
      do i = 2-mbc,mx+mbc
         i1 = i-2+imp                   !#  =  i-1 for amdq,  i for apdq

         do mq = 1,meqn
            sm = aux2(i1,1 + iface)
            sp = aux3(i1,1 + iface)

            bmasdq(i,mq) = min(sm,0.d0)*asdq(i,mq)
            bpasdq(i,mq) = max(sp,0.d0)*asdq(i,mq)
         enddo
      enddo

      return
      end
