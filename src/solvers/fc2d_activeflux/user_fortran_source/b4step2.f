      subroutine activeflux_b4step2(maxmx,maxmy,mbc,mx,my,
     &      meqn,q,xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer maxmx, maxmy, mbc, mx, my, meqn, maux
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      return
      end
