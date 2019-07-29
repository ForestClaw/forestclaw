      subroutine cudaclaw5_set_capacity(mx,my,mbc,dx,dy,
     &      area,mcapa,maux,aux)
      implicit none

      integer mbc, mx, my, maux, mcapa
      double precision dx, dy
      !!! Change arguments order
      double precision  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision dxdy

      dxdy = dx*dy

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            !!! Change arguments order
            aux(mcapa,i,j) = area(i,j)/dxdy
         enddo
      enddo

      end
