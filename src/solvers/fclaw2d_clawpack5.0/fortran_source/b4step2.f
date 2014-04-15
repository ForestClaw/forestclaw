c     ============================================
      subroutine b4step2(mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called from claw2 before each call to step2.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 
c
c     
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
c
      return
      end

