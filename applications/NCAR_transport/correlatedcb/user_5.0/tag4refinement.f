      subroutine clawpack5_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision qmin, qmax, dq,dq1,dq2

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1)
      qmax = q(mq,1,1)
      do j = 1,my
         do i = 1,mx
            dq1 = abs(q(mq,i+1,j) - q(mq,i-1,j))
            dq2 = abs(q(mq,i,j+1) - q(mq,i,j-1))
            dq = max(dq1,dq2)
            if (dq .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
