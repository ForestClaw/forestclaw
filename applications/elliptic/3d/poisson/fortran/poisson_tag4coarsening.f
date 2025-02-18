subroutine tag4coarsening(mx,my,mz,mbc,meqn,
     &      xlower,ylower,zlower,dx,dy,dz, blockno, 
     &      q0, q1, q2, q3, q4, q5, q6, q7,
     &      coarsen_threshold, initflag, tag_patch)
      implicit none

      integer mx,my,mz, mbc, meqn, tag_patch, initflag
      integer blockno
      double precision xlower(0:3), ylower(0:3), zlower(0:3)
      double precision dx, dy, dz
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q4(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q5(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q6(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      double precision q7(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

      integer mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1
      mq = 1
      qmin = q0(1,1,1,mq)
      qmax = q0(1,1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q0,qmin,qmax, coarsen_threshold,initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &              mq,q1,qmin,qmax, coarsen_threshold,
     &              initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &              mq,q2,qmin,qmax,coarsen_threshold,
     &              initflag,tag_patch )
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q3,qmin,qmax,coarsen_threshold,
     &      initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q4,qmin,qmax,coarsen_threshold,
     &      initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q5,qmin,qmax,coarsen_threshold,
     &      initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q6,qmin,qmax,coarsen_threshold,
     &      initflag,tag_patch)
      if (tag_patch == 0) return

      call poisson_get_minmax(mx,my,mz,mbc,meqn,
     &      mq,q7,qmin,qmax,coarsen_threshold,
     &      initflag,tag_patch)

      end

      subroutine poisson_get_minmax(mx,my,mz,mbc,meqn,mq,q,
     &      qmin,qmax,coarsen_threshold,initflag,tag_patch)

      implicit none
      integer mx,my,mz,mbc,meqn,mq,tag_patch,initflag
      double precision coarsen_threshold
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
      integer i,j,k

      do k = 1-mbc,mz+mbc
         do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
               if (abs(q(i,j,k,1)) .gt. coarsen_threshold) then
c                 # We won't coarsen this family because at least one
c                 # grid fails the coarsening test.
                  tag_patch = 0
                  return
               endif
            enddo
         enddo
      enddo

      end
