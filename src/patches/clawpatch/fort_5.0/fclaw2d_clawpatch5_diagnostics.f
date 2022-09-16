c    # ----------------------------------------------------------------------------------
c>     @file
c>     Diagnostics
c    # ----------------------------------------------------------------------------------

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_conscheck_t
c>
c> Implementation for clawpack 5
c>
c> @details @copydetails ::clawpatch_fort_conscheck_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_conservation_check(mx,my,
     &      mbc,meqn, dx,dy,area,q,sum)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy
      double precision sum(meqn)
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

c      include 'fclaw2d_metric_terms.i'

      integer i,j,m
      double precision dxdy
      integer*8 cont, fclaw_map_get_context
      logical fclaw2d_map_is_used

      cont = fclaw_map_get_context()

      dxdy = dx*dy
      do m = 1,meqn
         if (fclaw2d_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  sum(m) = sum(m) + q(m,i,j)*area(i,j)
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  sum(m) = sum(m) + q(m,i,j)*dx*dy
               enddo
            enddo
         endif
      enddo

      end

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_area_t
c>
c> Implementation for clawpack 5
c>
c> @details @copydetails ::clawpatch_fort_area_t
c--------------------------------------------------------------------
      double precision function
     &      fclaw2d_clawpatch5_fort_compute_patch_area(mx,my,
     &      mbc,dx,dy,area)
      implicit none

      integer mx,my, mbc
      double precision dx, dy
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

c      include 'fclaw2d_metric_terms.i'

      integer i,j
      double precision sum
      integer*8 cont, fclaw_map_get_context
      logical fclaw2d_map_is_used

      cont = fclaw_map_get_context()

      if (fclaw2d_map_is_used(cont)) then
         sum = 0
         do j = 1,my
            do i = 1,mx
               sum = sum + area(i,j)
            enddo
         enddo
      else
         sum = dx*dy*mx*my
      endif

      fclaw2d_clawpatch5_fort_compute_patch_area = sum

      end

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_norm_t
c>
c> Implementation for clawpack 5
c>
c> @details @copydetails ::clawpatch_fort_norm_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_compute_error_norm(
     &      blockno,mx,my, mbc,meqn, dx,dy,area,error,error_norm)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy
      double precision error_norm(meqn,3)
      double precision error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

c      include 'fclaw2d_metric_terms.i'

      integer i,j,m
      double precision dxdy, eij
      integer*8 cont, fclaw_map_get_context
      logical fclaw2d_map_is_used

      cont = fclaw_map_get_context()

c     # error_norm(:) comes in with values;  do not initialize it here!
      dxdy = dx*dy
      do m = 1,meqn
         if (fclaw2d_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  eij = abs(error(m,i,j))
                  error_norm(m,1) = error_norm(m,1) +
     &                  eij*area(i,j)
                  error_norm(m,2) = error_norm(m,2) +
     &                  eij**2*area(i,j)
                  error_norm(m,3) = max(eij,error_norm(m,3))
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  eij = abs(error(m,i,j))
                  error_norm(m,1) = error_norm(m,1) +
     &                  eij*dxdy
                  error_norm(m,2) = error_norm(m,2) +
     &                  eij**2*dxdy
                  error_norm(m,3) = max(eij,error_norm(m,3))
               enddo
            enddo
         endif
      enddo

      end
