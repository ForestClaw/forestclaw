c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
      subroutine compute_sum(mx,my,mbc,meqn,dx,dy,q,area,sum)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, sum
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      include 'metric_terms.i'

      integer i,j
      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         sum = 0
         do i = 1,mx
            do j = 1,my
               sum = sum + q(i,j,1)*area(i,j)
            enddo
         enddo
      else
         sum = 0
         do i = 1,mx
            do j = 1,my
               sum = sum + q(i,j,1)*dx*dy
            enddo
         enddo
      endif
      end
