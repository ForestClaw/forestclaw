      subroutine tag4coarsening_dq(mx,my,mbc,meqn,
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
         call get_minmax(mx,my,mbc,meqn,mq,q0,dq,
     &         coarsen_threshold,tag_patch)
         if (tag_patch .eq. 0) return

         call get_minmax(mx,my,mbc,meqn,mq,q1,dq,
     &         coarsen_threshold,tag_patch)
         if (tag_patch .eq. 0) return

         call get_minmax(mx,my,mbc,meqn,mq,q2,dq,
     &         coarsen_threshold,tag_patch)
         if (tag_patch .eq. 0) return

         call get_minmax(mx,my,mbc,meqn,mq,q3,dq,
     &         coarsen_threshold,tag_patch)
         if (tag_patch .eq. 0) return
      enddo

c      if (dq .gt. coarsen_threshold) then
c         tag_patch = 1
c      endif


      end

      subroutine get_minmax(mx,my,mbc,meqn,mq,q,dq,
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold
      double precision qmin,qmax,dq
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j
      double precision dqi, dqj

      do i = 1,mx
         do j = 1,my
            dqi = abs(q(i+1,j,mq) - q(i-1,j,mq))
            dqj = abs(q(i,j+1,mq) - q(i,j-1,mq))
            dq  = max1(dq,dqi, dqj)
            if (dq .gt. coarsen_threshold) then
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
