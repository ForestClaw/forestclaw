      subroutine heat_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, initflag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, initflag
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
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call heat_get_minmax(blockno, mx,my,mbc,meqn,
     &      mq,q0,dx,dy, coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(blockno, mx,my,mbc,meqn,
     &              mq,q1,dx,dy, coarsen_threshold,initflag,
     &              tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(blockno, mx,my,mbc,meqn,
     &              mq,q2,dx,dy,coarsen_threshold,initflag,
     &              tag_patch)
      if (tag_patch == 0) return

      call heat_get_minmax(blockno, mx,my,mbc,meqn,
     &      mq,q3,dx,dy,coarsen_threshold,initflag,
     &      tag_patch)

      end

      subroutine heat_get_minmax(blockno, mx,my,mbc,meqn,mq,q,
     &      dx,dy,coarsen_threshold,initflag,tag_patch)
      implicit none

      integer mx,my,mbc,meqn,mq,tag_patch,initflag
      integer blockno
      double precision coarsen_threshold, dx,dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision qvec(3), dmax, heat_eval_refinement

      integer i,j


      if (initflag .ne. 0) then
         tag_patch = 0
         return
      endif

      do i = 1,mx+1
         do j = 1,my+1
            qvec(1) = q(i,j,mq)
            qvec(2) = q(i-1,j,mq)
            qvec(3) = q(i,j-1,mq)

            dmax = heat_eval_refinement(qvec,dx,dy)

            if (abs(dmax) .gt. coarsen_threshold) then
               tag_patch = 0
               return
            endif

         enddo
      enddo

      end
