      subroutine initialize(mx,my,meqn,mbc,
     &      xlower,ylower,dx,dy,q,curvature,area)
     &      bind(c,name="initialize")
      implicit none

      integer meqn, mbc, mx, my, blockno
      double precision xlower, ylower, dx, dy, dxdy
      integer*8 cont
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      integer i, j

      include 'fclaw2d_metric_terms.i'

      dxdy = dx*dy
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            q(i,j,1) = area(i,j)/dxdy
            q(i,j,2) = curvature(i,j)
         enddo
      enddo

      return
      end
