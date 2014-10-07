      subroutine b4step2_manifold(mbc,mx,my,
     &      dx,dy,time,maux,aux,blockno,xd,yd,zd)
      implicit none

      integer mbc, mx, my, maux, blockno
      double precision xlower, ylower, dx, dy, time
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      call compute_velocity_psi(mx,my,mbc,dx,dy,
     &      time,blockno,xd,yd,zd,aux,maux)


      return
      end
