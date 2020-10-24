      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax, qij, qipj, qimj, qijp,qijm
      double precision qlap, dx2, dy2

      tag_patch = 0
      dx2 = dx*dx
      dy2 = dy*dy

c     # Refine based only on first variable in system.
      mq = 1
      do j = 1,my
         do i = 1,mx
            qij = q(i,j,mq)
            qipj = q(i+1,j,mq)
            qimj = q(i-1,j,mq)
            qijp = q(i,j+1,mq)
            qijm = q(i,j-1,mq)
            qlap = (qipj - 2*qij + qimj)/dx2 + 
     &             (qijp - 2*qij + qijm)/dy2
            if (abs(qlap) .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
