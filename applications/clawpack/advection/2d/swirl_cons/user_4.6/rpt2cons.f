      subroutine rpt2cons(ixy,maxm,meqn,mwaves,mbc,mx,
     &                       ql,qr,aux1,aux2,aux3,imp,asdq,
     &                       bmasdq,bpasdq)
      
      implicit none

      integer ixy, icoor, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)


      integer iuv, iface,i,j,i1
      double precision vl,vr,vhat

      iface = 3-ixy
      do i = 2-mbc, mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

c        # Lower face (at imp faces)
         vl = aux1(i1,iface)
         vr = aux2(i1,iface)
         vhat = (vl+vr)/2.d0

         bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)

c        # Upper face (at imp faces)
         vl = aux2(i1,iface)
         vr = aux3(i1,iface)
         vhat = (vl + vr)/2.d0
         bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)
      enddo


      return
      end
