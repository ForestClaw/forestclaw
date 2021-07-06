      subroutine clawpack46_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
  
      if (init_flag .eq. 1) then
          tag_patch = 0
          return
      endif


      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user46_get_minmax(blockno, mx,my,mbc,meqn,mq,q0, 
     &                       qmin,qmax, dx,dy,xlower(0),ylower(0),
     &                       coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(blockno, mx,my,mbc,meqn,mq,q1,
     &                       qmin,qmax, dx,dy,xlower(1),ylower(1),
     &                       coarsen_threshold, tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(blockno, mx,my,mbc,meqn,mq,q2,
     &                       qmin,qmax, dx,dy,xlower(2),ylower(2),
     &                       coarsen_threshold, tag_patch)
      if (tag_patch == 0) return

      call user46_get_minmax(blockno, mx,my,mbc,meqn,mq,q3,
     &                       qmin,qmax, dx,dy, xlower(3),ylower(3),
     &                       coarsen_threshold, tag_patch)

      end

      subroutine user46_get_minmax(blockno, mx,my,mbc,meqn,mq,q,
     &      qmin,qmax, dx,dy,xlower, ylower, 
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, blockno
      double precision coarsen_threshold, dx, dy, xlower,ylower
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision pi
      common /compi/ pi

      integer i,j, ii,jj
      double precision quad(-1:1,-1:1), xc, yc

      logical exceeds_th, minmax_exceeds_th

      do i = 1,mx
         do j = 1,my
             xc = xlower + (i-0.5)*dx
             yc = ylower + (j-0.5)*dy
             qmin = min(qmin,q(i,j,mq))
             qmax = max(qmax,q(i,j,mq))
              do ii = -1,1
                  do jj = -1,1
                        quad(ii,jj) = q(i+ii,j+jj,mq)
                  end do
              end do
              exceeds_th = minmax_exceeds_th(blockno,
     &                     q(i,j,mq),qmin,qmax,quad, dx,dy,xc,yc, 
     &                     coarsen_threshold)

            if (exceeds_th) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
