      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, initflag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, initflag
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq
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

      call heat_get_minmax(mx,my,mbc,meqn,
     &      mq,q0,qmin,qmax, coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(mx,my,mbc,meqn,
     &              mq,q1,qmin,qmax, coarsen_threshold,
     &              initflag,tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(mx,my,mbc,meqn,
     &              mq,q2,qmin,qmax,coarsen_threshold,initflag,
     &              tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(mx,my,mbc,meqn,
     &      mq,q3,qmin,qmax,coarsen_threshold,initflag,tag_patch)

      end

      subroutine heat_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,initflag,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, initflag
      double precision coarsen_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      
      integer i,j
      double precision qmin, qmax, qlow, qhi

      qlow = -1
      qhi = 1

      do i = 1,mx
         do j = 1,my
            if (initflag .ne. 0) then
                qmin = min(q(i,j,1),qmin)
                qmax = max(q(i,j,1),qmax)
                if (qmax-qmin .gt. coarsen_threshold) then
                   tag_patch = 0
                   return
                endif            
            else
                if (q(i,j,1) .gt. qlow + coarsen_threshold .and. 
     &               q(i,j,1) .lt. qhi - coarsen_threshold) then
c                   # We won't coarsen this family because at least one
c                   # grid fails the coarsening test.
                    tag_patch = 0
                    return
                endif
            endif
         enddo
      enddo

      end
