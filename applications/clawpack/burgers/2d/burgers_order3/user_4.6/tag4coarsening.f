c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine tag4coarsening(mx,my,mbc,meqn,
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

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1

      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call periodic_get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      dx,dy,xlower(0),ylower(0),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      dx,dy,xlower(1),ylower(1),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      dx,dy,xlower(2),ylower(2),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      dx,dy,xlower(3),ylower(3),coarsen_threshold,tag_patch)

      end

      subroutine periodic_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,dx,dy,xlower,ylower, 
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch
      double precision coarsen_threshold
      double precision dx,dy,xlower,ylower
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer example
      common /example_comm/ example  


      double precision qmin,qmax
      double precision qx, qy, xc,yc
      integer i,j

      do i = 1,mx
         do j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            if (example .eq. 0) then
               qmin = min(q(i,j,mq),qmin)
               qmax = max(q(i,j,mq),qmax)
               if (qmax-qmin .gt. coarsen_threshold) then

                  tag_patch = 0
                  return
               endif            
c               qx = (q(i+1,j,1)-q(i-1,j,1))/(2*dx)
c               qy = (q(i,j+1,1)-q(i,j-1,1))/(2*dy)
c               if (abs(qx) .gt. coarsen_threshold .or.
c     &               abs(qy) .gt. coarsen_threshold) then
c                  tag_patch = 0
c                  return
c               endif
            elseif (example .eq. 1) then
               if (abs(q(i,j,mq)) .gt. coarsen_threshold) then
                   tag_patch = 0
                   return
               endif
            elseif (example .eq. 2) then
               qx = (q(i+1,j,1)-q(i-1,j,1))/(2*dx)
               qy = (q(i,j+1,1)-q(i,j-1,1))/(2*dy)
               if (abs(qx) .gt. coarsen_threshold .or.
     &               abs(qy) .gt. coarsen_threshold) then
                  tag_patch = 0
                  return
               endif
            endif

         enddo
      enddo

      end
