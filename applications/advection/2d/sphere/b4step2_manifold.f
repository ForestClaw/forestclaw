      subroutine b4step2_manifold(maxmx,maxmy,mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux,
     &      xp,yp,zp,xd,yd,zd)
      implicit none

      integer maxmx, maxmy, mbc, mx, my, meqn, maux
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      call compute_velocity_psi(mx,my,mbc,dx,dy,
     &      time,xd,yd,zd,aux,maux)


      return
      end
