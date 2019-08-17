      subroutine torus_tag4coarsening(mx,my,mbc,meqn,
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

      tag_patch = 0
      return

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
     &      mq,q0,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(1), ylower(1), dx,dy,
     &      mq, q1,qmin,qmax,  coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(2), ylower(2), dx,dy,
     &      mq,q2,qmin,qmax, coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user_get_minmax(blockno, mx,my,mbc,meqn,
     &      xlower(3), ylower(3), dx,dy,
     &      mq,q3,qmin,qmax, coarsen_threshold,tag_patch)

      end

      subroutine user_get_minmax(blockno, mx,my,mbc,meqn, 
     &      xlower, ylower, dx,dy, mq,q,
     &      qmin,qmax,tag_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, blockno
      double precision tag_threshold, t, xlower, ylower, dx,dy
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer example
      common /example_comm/ example

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer*8 cont, get_context

      double precision xc,yc, r, xp, yp,zp
      logical refine
      integer i,j

c     # In case we need to call the physical mapping
      cont = get_context()

      refine = .false.
      do i = 1,mx
          do j = 1,my

              if (refine_pattern .eq. 0) then
                  refine = q(i,j,mq) .gt.  tag_threshold              
              else
                  xc = xlower + (i-0.5)*dx
                  yc = ylower + (j-0.5)*dy
                  call fclaw2d_map_c2m(cont,blockno,xc,yc,
     &                                     xp, yp,zp);
                  if (refine_pattern .eq. 1) then
                      refine = xp .lt. 0
                  elseif (refine_pattern .eq. 2) then
                      r = sqrt(xp**2 + yp**2)
                      refine = r .gt. 1.d0
                  endif
              endif
              if (refine) then
                 tag_patch = 0
                 return
              endif
          enddo
      enddo

      end
