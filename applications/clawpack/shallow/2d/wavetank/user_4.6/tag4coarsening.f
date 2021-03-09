c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine wavetank_fort_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, time, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy, time
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

      call user_get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      xlower(0),ylower(0),dx,dy,time, 
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      xlower(1),ylower(1),dx,dy, time, 
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      xlower(2),ylower(2),dx,dy, time, 
     &      coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      xlower(3),ylower(3),dx,dy, time, 
     &      coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax, xlower,ylower,dx,dy, time, 
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision xlower, ylower, dx,dy, time
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision b, bathy, xc, hij
      integer i,j

      tag_patch = 0

      do i = 1,mx
          xc = xlower + (i-0.5)*dx          
          if (time .ge. 15 .and. xc .lt.  -140.0) then
              tag_patch = 1
              return
          else if (time .ge. 30 .and. xc .lt. -100) then
              tag_patch = 1
              return
          elseif (time .ge. 45 .and. xc .lt. -70) then
              tag_patch = 1
              return
          elseif (time .ge. 60 .and. xc .lt. -50) then
              tag_patch = 1
             return
          endif
          b = bathy(xc,0)
          hij = q(i,j,1)
          qmin = min(hij + b,qmin)
          qmax = max(hij + b,qmax)
          if (hij .gt. 0 .and. abs(qmax) .gt. coarsen_threshold) then
c             # We won't coarsen this family because at least one
c             # grid fails the coarsening test.
              tag_patch = 0
              return
          endif
      enddo

      end
