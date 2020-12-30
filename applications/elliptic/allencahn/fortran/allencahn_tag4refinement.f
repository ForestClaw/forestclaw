      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax, qlow, qhi

      tag_patch = 0

c     # Refine based only on first variable in system.
      qlow = -1
      qhi = 1
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (init_flag .ne. 0) then
               if (qmax-qmin .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            elseif (q(i,j,1) .gt. qlow + tag_threshold .and. 
     &          q(i,j,1) .lt. qhi - tag_threshold) then
                tag_patch = 1
                return
            endif
         enddo
      enddo

      end
