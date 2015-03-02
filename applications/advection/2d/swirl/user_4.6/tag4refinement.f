      subroutine tag4refinement(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,blockno,q,
     &      refine_threshold, init_flag,
     &      tag_for_refinement)
      implicit none

      integer mx,my, mbc, meqn, tag_for_refinement, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision refine_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision xc,yc, qmin, qmax

      tag_for_refinement = 0

c     # Refine based only on first variable in system.
      do mq = 1,meqn
         qmin = 100.d0
         qmax = -100.d0
         do i = 1,mx
            do j = 1,my
               if (init_flag .eq. 1) then
                  xc = xlower + (i-0.5)*dx
                  yc = ylower + (j-0.5)*dy
                  if (abs(xc-0.5) .lt. dy) then
                     tag_for_refinement = 1
                     return
                  endif
               else
c                 # Exit immediately if the refinement criteria is met
                  qmin = min(q(i,j,mq),qmin)
                  qmax = max(q(i,j,mq),qmax)
                  if (qmax - qmin .gt. refine_threshold) then
                     tag_for_refinement = 1
                     return
                  endif
               endif
            enddo
         enddo
      enddo

      end

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_for_coarsening)
      implicit none

      integer mx,my, mbc, meqn, tag_for_coarsening
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      tag_for_coarsening = 0
      qmin = 100.d0
      qmax = -100.d0
      call tag_get_minmax(mx,my,mbc,meqn,q0,qmin,qmax)
      call tag_get_minmax(mx,my,mbc,meqn,q1,qmin,qmax)
      call tag_get_minmax(mx,my,mbc,meqn,q2,qmin,qmax)
      call tag_get_minmax(mx,my,mbc,meqn,q3,qmin,qmax)
      if (qmax - qmin .lt. coarsen_threshold) then
         tag_for_coarsening = 1
         return
      endif

      end

      subroutine tag_get_minmax(mx,my,mbc,meqn,q,
     &      qmin,qmax)

      implicit none
      integer mx,my,mbc,meqn
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j,mq

      mq = 1
      do i = 1,mx
         do j = 1,my
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
         enddo
      enddo

      end
