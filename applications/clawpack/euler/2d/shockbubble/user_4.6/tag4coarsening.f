      subroutine clawpack46_tag4coarsening(mx,my,mbc,meqn,
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
      double precision qmin, qmax, dq

      tag_patch = 0
      dq = 0
      do mq = 1,1
         call clawpack46_get_minmax(mx,my,mbc,meqn,mq,q0,dq)
         call clawpack46_get_minmax(mx,my,mbc,meqn,mq,q1,dq)
         call clawpack46_get_minmax(mx,my,mbc,meqn,mq,q2,dq)
         call clawpack46_get_minmax(mx,my,mbc,meqn,mq,q3,dq)
      enddo
      if (dq .lt. coarsen_threshold) then
         tag_patch = 1
         return
      endif

      end

      subroutine clawpack46_get_minmax(mx,my,mbc,meqn,mq,q,dq)

      implicit none
      integer mx,my,mbc,meqn,mq
      double precision qmin,qmax,dq
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j
      double precision dqi, dqj

      do i = 1,mx
         do j = 1,my
c           qmin = min(q(i,j,mq),qmin)
c           qmax = max(q(i,j,mq),qmax)
            dqi = dabs(q(i+1,j,mq) - q(i-1,j,mq))
            dqj = dabs(q(i,j+1,mq) - q(i,j-1,mq))
            dq  = dmax1(dq,dqi, dqj)
         enddo
      enddo

      end
