      subroutine sphere_compute_error_fort(mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t, q,error,xp,yp,zp)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j
      double precision gaussian_sum

      include 'fclaw2d_metric_terms.i'

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
            error(i,j,1) = q(i,j,1) -
     &            gaussian_sum(xp(i,j),yp(i,j),zp(i,j))
         enddo
      enddo
      end
