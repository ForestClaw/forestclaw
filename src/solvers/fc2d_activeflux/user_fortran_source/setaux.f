      subroutine activeflux_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer maxmx, maxmy,mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)


      return
      end
