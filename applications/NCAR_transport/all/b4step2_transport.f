      subroutine b4step2_transport(mx,my,mbc,
     &      dx,dy,time,maux,aux,xd,yd,zd)
      implicit none

      integer mbc, mx, my, maux
      double precision xlower, ylower, dx, dy, time
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      call compute_velocity_psi(mx,my,mbc,dx,dy,
     &      time,xd,yd,zd,aux,maux)


      return
      end
