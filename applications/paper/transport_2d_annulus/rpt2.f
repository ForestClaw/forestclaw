      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)


c     # Must use edge velocities
      integer color_equation
      common /eqn_comm/ color_equation

      if (color_equation .eq. 1) then
         call rpt2adv_manifold(ixy,imp,maxm,meqn,mwaves,maux,
     &               mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      else
         call rpt2cons_manifold(ixy,imp,maxm,meqn,mwaves,maux,
     &               mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      endif


      end


