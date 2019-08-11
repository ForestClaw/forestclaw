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
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      c0, c1, c2, c3,
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

      double precision   c0(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   c1(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   c2(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   c3(-mbc:mx+mbc+1,-mbc:my+mbc+1)


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

      call user_get_minmax(mx,my,mbc,meqn,mq,q0,c0,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q1,c1,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q2,c2,qmin,qmax,
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q3,c3,qmin,qmax,
     &      coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(mx,my,mbc,meqn,mq,q,curvature,
     &      qmin,qmax,tag_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision tag_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer example
      common /example_comm/ example

      logical refine
      double precision refine_curvature

      integer i,j

      refine_curvature = tag_threshold


      do i = 1,mx
         do j = 1,my
             if (abs(curvature(i,j)-1) .gt. 
     &           refine_curvature) then
                 tag_patch = 0
                 return
             endif
         enddo
      enddo

      end
