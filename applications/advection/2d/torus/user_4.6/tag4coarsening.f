c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
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
      double precision qmin, qmax

      qmin = 0.025d0
      qmax = 0.975d0
      tag_patch = 1
      call torus46_tag_get_minmax(mx,my,mbc,meqn,q0,dx,dy,
     &      qmin,qmax,coarsen_threshold,tag_patch)
      if (tag_patch .eq. 0) then
         return
      endif
      call torus46_tag_get_minmax(mx,my,mbc,meqn,q1,dx,dy,
     &      qmin,qmax,coarsen_threshold,tag_patch)
      if (tag_patch .eq. 0) then
         return
      endif
      call torus46_tag_get_minmax(mx,my,mbc,meqn,q2,dx,dy,
     &      qmin,qmax,coarsen_threshold,tag_patch)
      if (tag_patch .eq. 0) then
         return
      endif
      call torus46_tag_get_minmax(mx,my,mbc,meqn,q3,dx,dy,
     &      qmin,qmax,coarsen_threshold,tag_patch)
      if (tag_patch .eq. 0) then
         return
      endif
      return

      end

      subroutine torus46_tag_get_minmax(mx,my,mbc,meqn,q,
     &      dx,dy,qmin,qmax,coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn, tag_patch
      double precision dx,dy,qx,qy,coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j,mq

      mq = 1
      do i = 1,mx
         do j = 1,my
            qx = (q(i+1,j,1)-q(i-1,j,1))/(2*dx)
            qy = (q(i,j+1,1)-q(i,j-1,1))/(2*dy)
            if (abs(qx) .gt. coarsen_threshold .or.
     &            abs(qy) .gt. coarsen_threshold) then
               tag_patch = 0
               return
            endif
c            if (q(i,j,1) .gt. qmin .and. q(i,j,1) .lt. 0.9) then
c               tag_patch = 0
c               return
c            endif
         enddo
      enddo

      end
