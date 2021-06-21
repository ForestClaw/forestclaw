      subroutine fclaw2d_clawpatch46_fort_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, initflag, tag_patch)
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

c     # Don't coarsen when initializing the mesh
      if (initflag .ne. 0) then
           tag_patch = 0
           return
      endif

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

      call fclaw2d_clawpatch46_get_minmax(mx,my,mbc,meqn,
     &      mq,q0,qmin,qmax, coarsen_threshold,initflag,
     &      tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_get_minmax(mx,my,mbc,meqn,
     &              mq,q1,qmin,qmax, coarsen_threshold,initflag,
     &              tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_get_minmax(mx,my,mbc,meqn,
     &              mq,q2,qmin,qmax,coarsen_threshold,initflag,
     &              tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_get_minmax(mx,my,mbc,meqn,
     &      mq,q3,qmin,qmax,coarsen_threshold,initflag,
     &      tag_patch)

      end

      subroutine fclaw2d_clawpatch46_get_minmax(mx,my,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,initflag,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, initflag
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      logical exceeds_th, fclaw2d_clawpatch46_exceeds_th

      integer i,j

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            exceeds_th = fclaw2d_clawpatch46_exceeds_th(
     &             q(i,j,mq),qmin,qmax,coarsen_threshold)
            if (exceeds_th) then
               tag_patch = 0
               return
            endif
c            if (qmax-qmin .gt. coarsen_threshold) then
cc              # We won't coarsen this family because at least one
cc              # grid fails the coarsening test.
c               tag_patch = 0
c               return
c            endif
         enddo
      enddo

      end
