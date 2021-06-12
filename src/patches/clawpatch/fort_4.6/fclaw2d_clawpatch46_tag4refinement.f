      subroutine fclaw2d_clawpatch46_fort_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      logical exceeds_th, fclaw2d_clawpatch46_exceeds_th

c     # Assume that we won't refine      
      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            exceeds_th = fclaw2d_clawpatch46_exceeds_th(
     &             q(i,j,mq),qmin,qmax,tag_threshold)
            if (exceeds_th) then
c              # Refine this patch               
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end

c     # check to see if value exceeds threshold
      logical function fclaw2d_clawpatch46_exceeds_th(
     &                 qval,qmin,qmax,threshhold)

      implicit none
      double precision qval,qmin,qmax,threshhold
      logical refine

      refine = .false.
      if (qval .gt. threshhold) then
         refine = .true.
      endif

      fclaw2d_clawpatch46_exceeds_th = refine

      end
