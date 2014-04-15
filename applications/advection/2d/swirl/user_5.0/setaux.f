c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
c     ============================================

c     # set auxiliary arrays
c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)


      implicit none

      integer mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i, j
      double precision xll, yll, psi

c     # constant velocities which are used if tperiod=0 is specified
c     # in setprob.data

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
         enddo
      enddo

      return
      end
