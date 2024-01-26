      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer maxmx, maxmy,mbc,mx,my,maux
      double precision xlower,ylower,dx,dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision yj

      ! do i = 1-mbc,mx+mbc
      !    do j = 1-mbc,my+mbc
      !       yj = ylower + (j-0.5d0)*dy
      !       aux(i,j,1) = yj
      !    enddo
      ! enddo

      return
      end
