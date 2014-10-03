      subroutine initialize(mx,my,meqn,mbc,cont,blockno,
     &      xlower,ylower,dx,dy,q,error,curvature,area)
     &      bind(c,name="initialize")
      implicit none

      integer meqn, mbc, mx, my, blockno
      double precision xlower, ylower, dx, dy, dxdy
      integer*8 cont
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision error(1-mbc:mx+mbc, 1-mbc:my+mbc)
      integer i, j
      logical fclaw2d_map_is_flat, isflat

      include 'metric_terms.i'

      dxdy = dx*dy
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            q(i,j,1) = area(i,j)
            q(i,j,2) = curvature(i,j)
            q(i,j,3) = error(i,j)
         enddo
      enddo

      return
      end
