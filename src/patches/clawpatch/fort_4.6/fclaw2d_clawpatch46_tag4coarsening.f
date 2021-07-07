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

      call fclaw2d_clawpatch46_test_refine(blockno,mx,my,mbc,meqn,
     &      mq,q0,qmin,qmax, dx,dy,xlower(0), ylower(0), 
     &      coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_test_refine(blockno,mx,my,mbc,meqn,
     &      mq,q1,qmin,qmax,dx,dy,xlower(1), ylower(1), 
     &      coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_test_refine(blockno,mx,my,mbc,meqn,
     &      mq,q2,qmin,qmax,dx,dy,xlower(2), ylower(2),
     &      coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch46_test_refine(blockno,mx,my,mbc,meqn,
     &      mq,q3,qmin,qmax,dx,dy,xlower(3), ylower(3),
     &      coarsen_threshold,initflag, tag_patch)

      end

      subroutine fclaw2d_clawpatch46_test_refine(blockno,mx,my,mbc,
     &      meqn,mq,q, qmin,qmax,dx,dy,xlower,ylower,
     &      coarsen_threshold,init_flag,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, init_flag, blockno
      double precision coarsen_threshold
      double precision qmin,qmax, dx, dy, xlower, ylower
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision xc,yc,quad(-1:1,-1:1),qval

      integer i,j, ii, jj

      logical exceeds_th, fclaw2d_clawpatch_exceeds_threshold
      logical(kind=4) :: is_ghost, fclaw2d_clawpatch46_is_ghost

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            qval = q(i,j,mq)
            is_ghost = fclaw2d_clawpatch46_is_ghost(i,j,mx,my)
            if (.not. is_ghost) then
               do ii = -1,1               
                  do jj = -1,1
                     quad(ii,jj) = q(i+ii,j+jj,mq)
                  end do
               end do
            endif
            exceeds_th = fclaw2d_clawpatch_exceeds_threshold(
     &             blockno, qval,qmin,qmax,quad, dx,dy,xc,yc,
     &             coarsen_threshold, init_flag, is_ghost)
            if (exceeds_th) then
c              # This patch exceeds coarsen threshold and so 
c              # should not be coarsened.   
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
