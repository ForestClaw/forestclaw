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
      double precision qmin, qmax, qlow, qhi, phi

      tag_patch = 0

c     # Refine based only on first variable in system.
      qlow = -1
      qhi = 1
      mq = 2
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
            phi = q(i,j,2)
            qmin = min(phi,qmin)
            qmax = max(phi,qmax)
            if (init_flag .ne. 0) then
               if (qmax-qmin .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            elseif (phi .gt. qlow + tag_threshold .and. 
     &              phi .lt. qhi - tag_threshold) then
                tag_patch = 1
                return
            endif
         enddo
      enddo

      end
