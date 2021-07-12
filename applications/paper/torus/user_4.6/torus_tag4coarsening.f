      subroutine torus_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      tag_patch = 0

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user_check_refine(blockno, mx,my,mbc,meqn,
     &      xlower(0), ylower(0), dx,dy,
     &      mq,q0,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_check_refine(blockno, mx,my,mbc,meqn,
     &      xlower(1), ylower(1), dx,dy,
     &      mq, q1,qmin,qmax,  coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_check_refine(blockno, mx,my,mbc,meqn,
     &      xlower(2), ylower(2), dx,dy,
     &      mq,q2,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_check_refine(blockno, mx,my,mbc,meqn,
     &      xlower(3), ylower(3), dx,dy,
     &      mq,q3,qmin,qmax, coarsen_threshold,tag_patch)

      end

      subroutine user_check_refine(blockno, mx,my,mbc,meqn, 
     &      xlower, ylower, dx,dy, mq,q,
     &      qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, blockno
      double precision coarsen_threshold, t, xlower, ylower, dx,dy
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      logical refine
      integer init_flag

      init_flag = 0
      refine = .false.

      call torus_tag4refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      blockno, q, coarsen_threshold, init_flag,refine)

      if (refine) then
          tag_patch = 0
      else
          tag_patch = 1
      endif

      end
