      subroutine clawpack46_tag4refinement(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,blockno, q,refine_threshold,
     &      init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision refine_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qmin, qmax
      double precision dq, dqi, dqj

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            if (abs(q(i,j,mq)) .gt. refine_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
