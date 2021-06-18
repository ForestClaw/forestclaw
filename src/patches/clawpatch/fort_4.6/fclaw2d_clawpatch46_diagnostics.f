c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
      subroutine fclaw2d_clawpatch46_fort_conservation_check
     &      (mx,my,mbc,mfields,dx,dy,area,q,sum,c_kahan)
      implicit none

      integer mx,my,mbc,mfields
      double precision dx, dy
      double precision sum(mfields), c_kahan(mfields)
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision t, y, dxdy

c      include 'fclaw2d_metric_terms.i'

      integer i,j,m
      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      cont = get_context()

      dxdy = dx*dy
      do m = 1,mfields
         if (fclaw2d_map_is_used(cont)) then
C            sum(m) = 0
c            c_kahan = 0
            do j = 1,my
               do i = 1,mx
                  y = q(i,j,m)*area(i,j) - c_kahan(m)
                  t = sum(m) + y
                  c_kahan(m) = (t-sum(m)) - y
                  sum(m) = t
c                  sum(m) = sum(m) + q(i,j,m)*area(i,j)
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  y = q(i,j,m)*dxdy - c_kahan(m)
                  t = sum(m) + y
                  c_kahan(m) = (t-sum(m)) - y
                  sum(m) = t
!!                  sum(m) = sum(m) + q(i,j,m)*dxdy
               enddo
            enddo
         endif
      enddo

      end

c     # Compute area of a patch
      double precision function
     &      fclaw2d_clawpatch46_fort_compute_patch_area
     &      (mx,my, mbc,dx,dy,area)
      implicit none

      integer mx,my, mbc
      double precision dx, dy
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

c      include 'fclaw2d_metric_terms.i'

      integer i,j
      integer*8 cont, get_context
      logical fclaw2d_map_is_used
      double precision sum

      cont = get_context()

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

      fclaw2d_clawpatch46_fort_compute_patch_area = sum

      end


      subroutine fclaw2d_clawpatch46_fort_compute_error_norm
     &   (blockno, mx,my,mbc,mfields,dx,dy,area,error,error_norm)
      implicit none

      integer mx,my,mbc,mfields, blockno
      double precision dx, dy
      double precision error_norm(mfields,3)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

c      include 'fclaw2d_metric_terms.i'

      integer i,j,m
      double precision dxdy, eij

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      cont = get_context()

c     # error_norm(:) comes in with values;  do not initialize it here!
      dxdy = dx*dy
      do m = 1,mfields
         if (fclaw2d_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  eij = abs(error(i,j,m))
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
                  eij = abs(error(i,j,m))
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
