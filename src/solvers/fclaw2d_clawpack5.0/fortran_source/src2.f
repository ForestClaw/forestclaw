c
c
c =========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &             q,maux,aux,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
c
c     # dummy source routine... does nothing
c
      return
      end
