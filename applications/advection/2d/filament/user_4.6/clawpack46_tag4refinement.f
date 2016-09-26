      subroutine tag4refinement_dq(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax
      double precision dq, dqi, dqj,xc,yc,rc

      tag_patch = 0

c     # Refine based only on first variable in system.
      qmin = q(1,1,1)
      qmax = q(1,1,1)
      dq = 0
      do j = 1,my
         do i = 1,mx
            if (init_flag .ne. 0) then
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               rc = sqrt((xc-0.5)**2 + (yc-1.0)**2)
               if ((0.25-2*dx) < rc .and. rc < (0.25 + 2*dx)) then
                  tag_patch = 1
                  return
               endif
            else

c              qmin = min(qmin,q(i,j,1))
c              qmax = max(qmax,q(i,j,1))
c              if (qmax-qmin .gt. tag_threshold) then
c                  tag_patch = 1
c                  return
c              endif

               dqi = dabs(q(i+1,j,1) - q(i-1,j,1))
               dqj = dabs(q(i,j+1,1) - q(i,j-1,1))
               dq  = dmax1(dq, dqi, dqj)
               if (dq .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            endif
         enddo
      enddo


      end
