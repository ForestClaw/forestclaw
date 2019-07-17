      subroutine rpt2cons_manifold(ixy,maxm,meqn,mwaves,mbc,mx,
     &                       ql,qr,aux1,aux2,aux3,imp,asdq,
     &                       bmasdq,bpasdq)
      
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)


      integer i, i1, k, idir, iface, m
      double precision gp, gm, vp, vm

      do i = 2-mbc, mx+mbc
          i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
          if (ixy .eq. 1) then
             gm = aux2(i1,7)      !! Bottom edge length
             gp = aux3(i1,7)      !! Top edge length
             vm = gm*aux2(i1,4)   !! cell-velocity projected to bottom edge
             vp = gp*aux2(i1,5)   !! cell-velocity projected to top edge
          else
             gm = aux2(i1,6)      !! Left edge length
             gp = aux3(i1,6)      !! Right edge length
             vm = gm*aux2(i1,2)   !! cell-velocity projected to left edge
             vp = gp*aux2(i1,3)   !! cell-velocity projected to right edge
          endif

          do m = 1,meqn
              bmasdq(i,m) = min(vm,0.d0)*asdq(i,m)
              bpasdq(i,m) = max(vp,0.d0)*asdq(i,m)
          end do
      enddo


      return
      end
