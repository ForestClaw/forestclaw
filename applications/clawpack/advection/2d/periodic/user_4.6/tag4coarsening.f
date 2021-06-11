c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine clawpack46_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.

      if (init_flag .ne. 0) then
            !! Don't coarsen while initializing the grids
            tag_patch = 0
            return
      endif

      tag_patch = 1   !! assume coarsening
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user_get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

      do i = 1,mx
         do j = 1,my
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if ((qmax-qmin) .gt. coarsen_threshold) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
