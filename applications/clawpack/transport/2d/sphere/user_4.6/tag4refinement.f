      subroutine sphere_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, curvature, tag_threshold, 
     &      init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision   curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer example
      common /example_comm/ example      

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      double precision curve_threshold


      logical refine_curve, refine


      cont = get_context()


      curve_threshold = tag_threshold
      tag_patch = 0

c     # Refine based only on first variable in system.
      refine = .false.
      refine_curve = .false.
      mq = 1
      do j = 1,my
         do i = 1,mx
            refine_curve = abs(curvature(i,j)-1.d0) 
     &                   .gt. curve_threshold
            if (initchoice .le. 2) then                
c                  refine = q(i,j,mq) .gt.  tag_threshold              
            else                
                write(6,'(A,A)') 'Refining not yet defined for ',
     &                  'example > 0'
                stop
            endif
            if (refine_curve) then
                tag_patch = 1
                return
            endif
         enddo
      enddo


      end

