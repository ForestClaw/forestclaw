c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================

c     # set auxiliary arrays
c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)


      implicit none

      integer maxmx, maxmy,mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision xll, yll, psi

c     # constant velocities which are used if tperiod=0 is specified
c     # in setprob.data

      return
      end
