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
      double precision qmin, qmax

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1)
      qmax = q(mq,1,1)
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            if (qmax - qmin .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
