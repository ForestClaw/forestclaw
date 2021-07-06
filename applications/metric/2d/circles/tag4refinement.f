      subroutine tag4refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      init_flag, blockno, tag_patch, xp, yp)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy, xc, yc, dr
      logical is_in_circle

      integer i,j, m

      include 'fclaw2d_metric_terms.i'

      tag_patch = 0
      do i = 1,mx
         do j = 1,my
            if (is_in_circle(xp(i,j),yp(i,j))) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine tag4coarsening(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      qcoarsened, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      double precision xlower, ylower, dx, dy
      double precision qcoarsened(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j
      double precision qmin, qmax

c     # The difference between this and the true "refinement" above is
c     # that we can't check ghost cells here.  Also, we may make the
c     # coarsening criteria different from the refinement criteria.
c     # Also, we don't check for an init_flag, since it is unlikely that
c     # we would coarsen an initial grid.

      qmin = 100.d0
      qmax = -100.d0
      tag_patch = 0
      do i = 1,mx
         do j = 1,my
            qmin = min(qcoarsened(i,j,1),qmin)
            qmax = max(qcoarsened(i,j,1),qmax)
            if (qmax - qmin .gt. 0.5d0) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
