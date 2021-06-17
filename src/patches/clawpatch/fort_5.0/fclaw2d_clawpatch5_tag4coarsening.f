c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine fclaw2d_clawpatch5_fort_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, initflag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch,initflag
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision q3(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(mq,1,1)
      qmax = q0(mq,1,1)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call fclaw2d_clawpatch5_get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,
     &      coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,
     &      coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,
     &      coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,
     &      coarsen_threshold,initflag,tag_patch)

      end

      subroutine fclaw2d_clawpatch5_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,initflag,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch,initflag
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer i,j

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            if (qmax - qmin .gt. coarsen_threshold) then
c              # We won't coarsen this family because at least one
c              # grid fails the coarsening test.
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
