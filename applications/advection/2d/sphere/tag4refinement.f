      subroutine tag_for_refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      q,init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qmin, qmax
      double precision dq, dqi, dqj

      qmin = 100.d0
      qmax = -100.d0
      tag_patch = 0
      do mq = 1,meqn
         do i = 1,mx
            do j = 1,my

               if (init_flag == 1 .and. .false.) then
c                 # Be careful : comp. coord. are [0,1]x[0,1]
                  xc = xlower + (i-0.5)*dx
                  yc = ylower + (j-0.5)*dy
                  if (abs(yc-0.5d0) < dx) then
                     tag_patch = 1
                     return
                  endif
c                  tag_patch = 0
c                  return
               else
c                  tag_patch = 1
c                  return
                  qmin = min(q(i,j,mq),qmin)
                  qmax = max(q(i,j,mq),qmax)
                  if (qmax - qmin .gt. 0.25d0) then
                     tag_patch = 1
                     return
                  endif
               endif
            enddo
         enddo
      enddo

      end
