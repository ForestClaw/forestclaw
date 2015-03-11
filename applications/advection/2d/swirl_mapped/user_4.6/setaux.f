c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &      dx,dy,maux,aux)
c     ============================================

c     # set auxiliary arrays
c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)


      implicit none

      integer  maxmx, maxmy, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i, j
      double precision xll, yll, psi, z0

c     # constant velocities which are used if tperiod=0 is specified
c     # in setprob.data

      z0 = 0.d0
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(i,j,1) = -(psi(xll, yll+dy,z0) - psi(xll,yll,z0)) / dy
            aux(i,j,2) =  (psi(xll+dx, yll,z0) - psi(xll,yll,z0)) / dx
         enddo
      enddo

      return
      end
