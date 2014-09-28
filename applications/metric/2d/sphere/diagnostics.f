      double precision function total_area(mx,my,mbc,area)
      implicit none

      integer mx,my,mbc
      double precision dx,dy
      double precision sum
      integer*8 cont, get_context
      integer blockno, get_block
      integer i,j

      include 'metric_terms.i'

      sum = 0
      do i = 1,mx
         do j = 1,my
            sum = sum + area(i,j)
         enddo
      enddo
      total_area = sum

      end


      subroutine min_grid_cell_area(mx,my,mbc,area,
     &      minvalue)
      implicit none

      integer mx,my,mbc
      double precision dx,dy
      double precision minvalue
      integer*8 cont, get_context
      integer blockno, get_block
      integer i,j

      include 'metric_terms.i'

c     # minvalue comes in with a value;  we only
c     # compare it here.
      do i = 1,mx
         do j = 1,my
            minvalue = min(minvalue,area(i,j))
         enddo
      enddo

      end

      subroutine max_grid_cell_area(mx,my,mbc,area,
     &      maxvalue)
      implicit none

      integer mx,my,mbc
      double precision dx,dy
      double precision maxvalue
      integer*8 cont, get_context
      integer blockno, get_block
      integer i,j

      include 'metric_terms.i'

c     # maxvalue comes in with a value;  we only
c     # compare it here.
      do i = 1,mx
         do j = 1,my
            maxvalue = max(maxvalue,area(i,j))
         enddo
      enddo

      end
