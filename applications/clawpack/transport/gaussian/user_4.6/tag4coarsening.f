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

c     # If we are initializing the refinement levels, don't coarsen.
      if (init_flag .eq. 1) then
          tag_patch = 0
          return
      endif

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user46_get_minmax(mx,my,mbc,meqn,mq,q0,dx,dy,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(mx,my,mbc,meqn,mq,q1,dx,dy,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(mx,my,mbc,meqn,mq,q2,dx,dy,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(mx,my,mbc,meqn,mq,q3,dx,dy,qmin,qmax,
     &      coarsen_threshold,tag_patch)

      end

      subroutine user46_get_minmax(mx,my,mbc,meqn,mq,q,
     &      dx,dy,qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold, dx,dy
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

      double precision pi
      common /compi/ pi

      double precision dq1,dq2,dq, dx2, dy2, s

      s = pi/2
      dx2 = 2*dx*s
      dy2 = 2*dy*s

      do i = 1,mx
         do j = 1,my
            dq1 = q(i+1,j,mq) - q(i-1,j,mq)/dx2
            dq2 = q(i,j+1,mq) - q(i,j-1,mq)/dy2
            dq = max(abs(dq1),abs(dq2))
            if (dq .gt. coarsen_threshold) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
