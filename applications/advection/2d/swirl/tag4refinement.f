      subroutine tag_for_refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      q,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qs(4), qmin, qmax

      qmin = 100.d0
      qmax = -100.d0
      tag_patch = 0
      do mq = 1,meqn
         do i = 1,mx
            do j = 1,my

               qmin = min(q(i,j,mq),qmin)
               qmax = max(q(i,j,mq),qmax)

               do m = 1,4
                  if (qs(m) > 0.5d0) then
c                    # We found one cell to refine
c                     tag_patch = 1
c                     return
                  endif
               enddo

c              xc = xlower + (i-0.5)*dx
c              yc = ylower + (j-0.5)*dy
c              if (abs(yc - 0.5d0) < dy) then
c                 tag_patch = 1
c              endif
            enddo
         enddo
      enddo
      if (qmax - qmin > 0.5d0) then
         tag_patch = 1
      endif

      end
