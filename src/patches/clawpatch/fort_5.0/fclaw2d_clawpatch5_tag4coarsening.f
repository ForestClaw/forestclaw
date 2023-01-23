c> @file
c> tag4coarsening routine for clawpack 5

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_tag4coarsening_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::clawpatch_fort_tag4coarsening_t
c--------------------------------------------------------------------
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
      double precision qmin(meqn), qmax(meqn)

c     # Don't coarsen when initializing the mesh
      if (initflag .ne. 0) then
           tag_patch = 0
           return
      endif

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      do mq = 1,meqn
         qmin(mq) = q0(mq,1,1)
         qmax(mq) = q0(mq,1,1)
      end do

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call fclaw2d_clawpatch5_test_refine(blockno,mx,my,mbc,meqn,
     &      q0,qmin,qmax, dx,dy,xlower(0), ylower(0), 
     &      coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_test_refine(blockno,mx,my,mbc,meqn,
     &      q1,qmin,qmax,dx,dy,xlower(1), ylower(1), 
     &      coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_test_refine(blockno,mx,my,mbc,meqn,
     &              q2,qmin,qmax,dx,dy,xlower(2), ylower(2),
     &              coarsen_threshold,initflag, tag_patch)
      if (tag_patch == 0) return

      call fclaw2d_clawpatch5_test_refine(blockno,mx,my,mbc,meqn,
     &      q3,qmin,qmax,dx,dy,xlower(3), ylower(3),
     &      coarsen_threshold,initflag, tag_patch)

      end

c--------------------------------------------------------------------
c> @brief Tests a single patch to see if it exceeds threshold
c>
c> @param[in] blockno the block number
c> @param[in] mx, my the number of cells in the x and y directions
c> @param[in] mbc the number of ghost cells
c> @param[in] meqn the number of equations
c> @param[in] mq the equation to test
c> @param[in] q the solution
c> @param[in,out] qmin the minimum value in q
c> @param[in,out] qmax the maximum value in q
c> @param[in] dx, dy the spacings in the x and y directions
c> @param[in] xlower, ylower the lower left coordinate of the patch
c> @param[in] coarsen_threshold the threshold
c> @param[in] init_flag true if in init stage
c> @param[in,out] tag_patch passed in as [1] may be set to [0] if it
c>                should not be coarsened
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_test_refine(blockno, mx,my,mbc,
     &                meqn,q, qmin,qmax,dx,dy,xlower,ylower,
     &                coarsen_threshold,init_flag,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch,init_flag,blockno
      double precision coarsen_threshold
      double precision qmin(mq),qmax(mq), dx, dy, xlower, ylower
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      double precision xc,yc,quad(-1:1,-1:1,meqn), qval(meqn)
      integer exceeds_th, fclaw2d_clawpatch_tag_criteria

      integer i,j, ii, jj

      logical(kind=4) is_ghost, fclaw2d_clawpatch5_is_ghost

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            do mq=1,meqn
               qval(mq) = q(mq,i,j)
               qmin(mq) = min(q(mq,i,j),qmin(mq))
               qmax(mq) = max(q(mq,i,j),qmax(mq))
            end do
            is_ghost = fclaw2d_clawpatch5_is_ghost(i,j,mx,my)
            if (.not. is_ghost) then
               do ii = -1,1
                  do jj = -1,1
                     do mq = 1,meqn
                        quad(ii,jj,mq) = q(mq,i+ii,j+jj)
                     end do                     
                  end do
               end do
            endif
            exceeds_th = fclaw2d_clawpatch_tag_criteria(
     &             blockno, qval,qmin,qmax,quad, dx,dy,xc,yc,
     &             coarsen_threshold, init_flag, is_ghost)

c           # -1 : Not conclusive (possibly ghost cell); don't tag for coarsening
c           # 0  : Does not pass threshold (tag for coarsening)      
c           # 1  : Passes threshold (do not tag for coarsening)
c           # exceeds_th = -1 is ignored here (logic of regridding seems
c           # to over-refine.  Not clear why.)
            if (exceeds_th .gt. 0) then
c              # This patch exceeds coarsen threshold and so 
c              # should not be coarsened.   
               tag_patch = 0
               return
            endif
         enddo
      enddo

      end
