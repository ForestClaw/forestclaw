      subroutine rpt2cons_cc(ixy,maxm,meqn,mwaves,mbc,mx,
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


      integer iuv, iface,i,j,i1,m, mq
      double precision vl,vr,v, nv

      iface = 3-ixy
      do i = 2-mbc, mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

         do mq = 1,meqn
c           # Get normal (times area) at this edge
c           # Cell centered velocity
            v = 0
c            do m = 1,3
cc              # Get normal for bmasdq
cc              # Normal at lower edge
c               nv = aux2(i1,1  + (iface-1)*3 + m)
c               vl = aux1(i1,10 + (mq-1)*3 + m)
c               vr = aux2(i1,10 + (mq-1)*3 + m)
c               v = v + 0.5*(vl + vr)*nv
c            enddo

            if (ixy .eq. 1) then
               vr = aux2(i1,4)
               vl = aux1(i1,5)
            else
               vr = aux2(i1,2)
               vl = aux1(i1,3)
            endif
            v = (vl+vr)/2.d0

            bmasdq(i,mq) = min(v,0.d0)*asdq(i,mq)

c           # Now do upper edge
c            v = 0
c            do m = 1,3
cc              # Get normal for bmasdq
cc              # Normal at lower edge
c               nv = aux3(i1,1+(iface-1)*3+m)
c               vl = aux2(i1,10+(mq-1)*3 + m)
c               vr = aux3(i1,10+(mq-1)*3 + m)
c               v = v + 0.5*(vl + vr)*nv
c            enddo

            if (ixy .eq. 1) then
               vr = aux3(i1,4)
               vl = aux2(i1,5)
            else
               vr = aux3(i1,2)
               vl = aux2(i1,3)
            endif
            v = (vl + vr)/2.d0
            bpasdq(i,mq) = max(v,0.d0)*asdq(i,mq)
         enddo
      enddo


      return
      end
