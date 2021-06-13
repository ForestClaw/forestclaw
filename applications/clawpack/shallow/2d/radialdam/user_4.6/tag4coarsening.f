c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine clawpatch46_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno, init_flag
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      if (init_flag .ne. 0) then
            tag_patch = 0
            return
      endif

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user_get_minmax(blockno, mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      xlower(0), ylower(0), dx, dy, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      xlower(1), ylower(1), dx, dy, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      xlower(2), ylower(2), dx, dy, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      xlower(3), ylower(3), dx, dy, coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(blockno, mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,xlower, ylower, dx, dy, 
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch,blockno
      double precision dx,dy,xlower, ylower
      double precision coarsen_threshold
      double precision quad(-1:1,-1:1)
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j, ii, jj

      logical exceeds_th, radialdam_exceeds_th
      double precision xc,yc

      do i = 2,mx-1         
         do j = 2,my-1  
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            do ii = -1,1
               do jj = -1,1
                  quad(ii,jj) = q(i+ii,j+jj,mq)
               end do
            end do

            exceeds_th = radialdam_exceeds_th(blockno, 
     &             q(i,j,mq),qmin,qmax,quad,dx,dy,xc,yc,
     &             coarsen_threshold)
            if (exceeds_th) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
