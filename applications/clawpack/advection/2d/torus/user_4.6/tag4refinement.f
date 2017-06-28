      subroutine clawpack46_tag4refinement(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,blockno,q,refine_threshold,
     &      init_flag,tag_for_refinement)
      implicit none

      integer mx,my, mbc, meqn, tag_for_refinement, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision refine_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision xc,yc, qmin, qmax,qx,qy

      tag_for_refinement = 0

c     # Refine based only on first variable in system.
      qmin = 0.025d0
      qmax = 0.975d0
      do mq = 1,meqn
         do i = 1,mx
            do j = 1,my
               qx = (q(i+1,j,1)-q(i-1,j,1))/(2*dx)
               qy = (q(i,j+1,1)-q(i,j-1,1))/(2*dy)
               if (abs(qx) .gt. refine_threshold .or.
     &               abs(qy) .gt. refine_threshold) then
                  tag_for_refinement = 1
                  return
               endif
c              # Tag grid if at least one value is in [qmin,qmax]
c               if (q(i,j,1) .gt. qmin
c     &               .and. q(i,j,1) .lt. qmax) then
c                  tag_for_refinement = 1
c                  return
c               endif
            enddo
         enddo
      enddo

      end
