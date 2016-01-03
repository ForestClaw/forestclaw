c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
      subroutine fclaw2d_fort_conservation_check(mx,my,mbc,meqn,
     &      dx,dy,area,q,sum)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, dxdy
      double precision sum(meqn)
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      include 'metric_terms.i'

      integer i,j,m
      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      cont = get_context()

      dxdy = dx*dy
      do m = 1,meqn
         if (fclaw2d_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  sum(m) = sum(m) + q(i,j,1)*area(i,j)
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  sum(m) = sum(m) + q(i,j,1)*dx*dy
               enddo
            enddo
         endif
      enddo

      end
