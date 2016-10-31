      subroutine clawpack46_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax
      double precision dq, dqi, dqj

      tag_patch = 0

c     # Refine based only on first variable in system.
      qmin = q(1,1,1)
      qmax = q(1,1,1)
      do j = 1,my
         do i = 1,mx
            dq = 0
            do mq = 1,1
                dqi = dabs(q(i+1,j,mq) - q(i-1,j,mq))
                dqj = dabs(q(i,j+1,mq) - q(i,j-1,mq))
                dq  = dmax1(dq, dqi, dqj)
                if (dq .gt. tag_threshold) then
                   tag_patch = 1
                   return
                endif
            enddo
         enddo
      enddo

      end
