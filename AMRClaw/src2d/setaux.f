c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c     # dummy routine when no auxiliary arrays
c
c
      implicit none
      integer maxmx, maxmy, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy

      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
c
       return
       end
