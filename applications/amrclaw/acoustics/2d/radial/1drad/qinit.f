c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c
      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc)
      dimension aux(maux, 1-mbc:mx+mbc)
c
      pi = 4.d0*datan(1.d0)
      width = 0.2d0
c
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         if (dabs(xcell-0.5d0) .le. width) then
             pressure = 1.d0 + dcos(pi*(xcell - 0.5d0)/width)
           else
             pressure = 0.d0
           endif
         q(1,i) = pressure
         q(2,i) = 0.d0
  150    continue
c
      return
      end
