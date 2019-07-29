      subroutine clawpatch_time_sync_setup(mx,my,mbc,dx,dy,
     &                                     area,edgelengths,
     &                                     area0,area1,area2,area3,
     &                                     el0, el1, el2, el3, 
     &                                     manifold)


      implicit none

      integer mx,my,mbc
      integer manifold
      double precision dx,dy

      double precision area0(my), area1(my), area2(mx), area3(mx)
      double precision el0(my), el1(my), el2(mx), el3(mx)

      integer i,j
      double precision dxdy

      include "metric_terms.i"

      if (manifold .eq. 1) then
         do j = 1,my
            area0(j) = area(1,j)
            area1(j) = area(mx,j)
            el0(j) = edgelengths(1,j,1)
            el1(j) = edgelengths(mx+1,j,1)
         enddo
    
         do i = 1,mx
            area2(i) = area(i,1)
            area3(i) = area(i,my)
            el2(i) = edgelengths(i,1,2)
            el3(i) = edgelengths(i,my+1,2)
         enddo
      else
         dxdy = dx*dy
         do j = 1,my
            area0(j) = dxdy
            area1(j) = dxdy
            el0(j) = dy
            el1(j) = dy
         enddo
         do i = 1,mx
            area2(i) = dxdy
            area3(i) = dxdy
            el2(i) = dx
            el3(i) = dx
         enddo
      endif


      end