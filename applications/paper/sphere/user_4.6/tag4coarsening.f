c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine sphere_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, t, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy, t
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(0), ylower(0), dx,dy,
     &      mq,q0,t, qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(1), ylower(1), dx,dy,
     &      mq, q1,t, qmin,qmax,  coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(2), ylower(2), dx,dy,
     &      mq,q2,t, qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(3), ylower(3), dx,dy,
     &      mq,q3,t, qmin,qmax, coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(blockno, mx,my,mbc,meqn, 
     &      xlower, ylower, dx,dy, mq,q,t, 
     &      qmin,qmax,tag_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, blockno
      double precision tag_threshold, t, xlower, ylower, dx,dy
      double precision qmin,qmax, lap
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer example
      common /example_comm/ example

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer*8 cont, get_context

      double precision xc,yc
      logical refine
      integer i,j

c     # In case we need to call the physical mapping
      cont = get_context()

      refine = .false.
      do i = 1,mx
         do j = 1,my
             xc = xlower + (i-0.5)*dx
             yc = ylower + (j-0.5)*dy
c            lap = (q(i+1,j,1) + q(i-1,j,1) + q(i,j+1,1) + q(i,j-1,1)
c     &            - 4*q(i,j,1))/dx**2
c            if (abs(lap) .gt. tag_threshold) then
c                  refine = .true.
c            endif
             if (initchoice .le. 3) then
                 refine = q(i,j,mq) .gt.  tag_threshold              
             endif
             if (refine) then
                 tag_patch = 0
                 return
             endif
         enddo
      enddo

      end
