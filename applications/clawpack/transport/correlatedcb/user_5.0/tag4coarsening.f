c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine clawpack5_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno, init_flag
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision qmin, qmax, dq, dq1, dq2

      if (init_flag .ne. 0) then
            tag_patch = 0
            return
      end if 

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(mq,1,1)
      qmax = q0(mq,1,1)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call user5_get_minmax(blockno, mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      dx,dy,xlower(0),ylower(0),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user5_get_minmax(blockno, mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      dx,dy,xlower(1),ylower(1),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user5_get_minmax(blockno, mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      dx,dy,xlower(2),ylower(2),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call user5_get_minmax(blockno, mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      dx,dy,xlower(3),ylower(3),coarsen_threshold,tag_patch)

      end

      subroutine user5_get_minmax(blockno, mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,dx,dy,xlower,ylower,coarsen_threshold,tag_patch)

      implicit none
      integer blockno, mx,my,mbc,meqn,mq,tag_patch
      double precision dx,dy,xlower, ylower
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer i,j

      double precision xc,yc,quad(-1:1,-1:1)
      integer ii,jj

      logical exceeds_th, transport_exceeds_th

      do i = 1,mx
          do j = 1,my
              xc = xlower + (i-0.5)*dx
              yc = ylower + (j-0.5)*dy
              qmin = min(qmin,q(mq,i,j))
              qmax = max(qmax,q(mq,i,j))
              do ii = -1,1
                  do jj = -1,1
                      quad(ii,jj) = q(mq,i+ii,j+jj)
                  end do
              end do
              exceeds_th = transport_exceeds_th(blockno,
     &                     q(mq,i,j),qmin,qmax,quad, dx,dy,xc,yc, 
     &                     coarsen_threshold)

              if (exceeds_th) then
c                 # We won't coarsen this family because at least one
c                 # grid fails the coarsening test.
                  tag_patch = 0
                  return
              endif
         enddo
      enddo

      end
