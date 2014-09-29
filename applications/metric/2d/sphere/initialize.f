      subroutine initialize(mx,my,meqn,mbc,cont,blockno,
     &      xlower,ylower,dx,dy,q,error,curvature,area)
     &      bind(c,name="initialize")
      implicit none

      integer meqn, mbc, mx, my, blockno
      double precision xlower, ylower, dx, dy
      integer*8 cont
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision error(1-mbc:mx+mbc, 1-mbc:my+mbc)
      integer i, j
      logical fclaw2d_map_is_flat, isflat

      include 'metric_terms.i'

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c            q(i,j,1) = area(i,j)
c            q(i,j,2) = curvature(i,j)
c            q(i,j,3) = error(i,j)
            q(i,j,1) = 0
            q(i,j,2) = 0
            q(i,j,3) = 0
         enddo
      enddo

      return
      end
