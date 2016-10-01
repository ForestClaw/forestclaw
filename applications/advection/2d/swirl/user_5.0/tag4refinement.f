      subroutine clawpack5_tag4refinement(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,blockno,q,
     &      refine_threshold, init_flag,
     &      tag_for_refinement)
      implicit none

      integer mx,my, mbc, meqn, tag_for_refinement, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision refine_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision xc,yc, qmin, qmax

      tag_for_refinement = 0

c     # Refine based only on first variable in system.
      do mq = 1,meqn
         qmin = 100.d0
         qmax = -100.d0
         do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
               if (init_flag .eq. 1) then
                  xc = xlower + (i-0.5)*dx
                  yc = ylower + (j-0.5)*dy
                  if (abs(xc-0.5) .lt. dy) then
                     tag_for_refinement = 1
                     return
                  endif
               else
c                 # Exit immediately if the refinement criteria is met
                  qmin = min(q(mq,i,j),qmin)
                  qmax = max(q(mq,i,j),qmax)
                  if (qmax - qmin .gt. refine_threshold) then
                     tag_for_refinement = 1
                     return
                  endif
               endif
            enddo
         enddo
      enddo

      end
