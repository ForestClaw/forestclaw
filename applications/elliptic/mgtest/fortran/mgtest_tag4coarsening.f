      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
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
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call mgtest_get_minmax(mx,my,mbc,meqn,
     &      mq,q0,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call mgtest_get_minmax(mx,my,mbc,meqn,
     &              mq,q1,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call mgtest_get_minmax(mx,my,mbc,meqn,
     &              mq,q2,qmin,qmax,coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call mgtest_get_minmax(mx,my,mbc,meqn,
     &      mq,q3,qmin,qmax,coarsen_threshold,tag_patch)

      end

      subroutine mgtest_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (qmax .gt. coarsen_threshold .or. 
     &              qmin .lt. -coarsen_threshold) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
