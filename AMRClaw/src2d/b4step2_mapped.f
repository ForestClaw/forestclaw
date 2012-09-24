      subroutine b4step2_mapped(mx,my,mbc,meqn,q,
     &      dx,dy,xp,yp,zp,xd,yd,zd,time,dt,maux,aux)
      implicit none

      integer maxmx, maxmy, mbc, mx, my, meqn, maux
      double precision dx,dy,time, dt
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      return
      end
