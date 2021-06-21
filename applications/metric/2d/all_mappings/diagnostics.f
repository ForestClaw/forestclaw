      double precision function total_area(mx,my,mbc,area)
      implicit none

      integer mx,my,mbc
      double precision dx,dy
      double precision sum
      integer i,j

      include 'fclaw2d_metric_terms.i'

      sum = 0
      do j = 1,my
         do i = 1,mx
            sum = sum + area(i,j)
         enddo
      enddo
      total_area = sum

      end


      subroutine min_grid_cell_area(mx,my,mbc,dx,dy,
     &      area,minvalue)
      implicit none

      integer mx,my,mbc
      double precision dx,dy,dxdy
      double precision minvalue
      integer i,j

      include 'fclaw2d_metric_terms.i'

c     # minvalue comes in with a value;  we only
c     # compare it here.
      dxdy = dx*dy
      do i = 1,mx
         do j = 1,my
            minvalue = min(minvalue,area(i,j)/dxdy)
         enddo
      enddo

      end

      subroutine max_grid_cell_area(mx,my,mbc,dx,dy,
     &      area,maxvalue)
      implicit none

      integer mx,my,mbc
      double precision dx,dy,dxdy
      double precision maxvalue
      integer i,j

      include 'fclaw2d_metric_terms.i'

c     # maxvalue comes in with a value;  we only
c     # compare it here.
      dxdy = dx*dy
      do i = 1,mx
         do j = 1,my
            maxvalue = max(maxvalue,area(i,j)/dxdy)
         enddo
      enddo

      end
