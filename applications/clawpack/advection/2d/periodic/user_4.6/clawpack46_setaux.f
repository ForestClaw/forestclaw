      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      double precision uvel, vvel
      common /comvelocity/ uvel, vvel

      integer i, j

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            
c           # difference stream function psi to get normal velocities:
            aux(i,j,1) = uvel
            aux(i,j,2) = vvel
         enddo
      enddo

      return
      end
