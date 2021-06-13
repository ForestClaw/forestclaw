      subroutine clawpack5_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq
      double precision qmin, qmax, quad(-1:1,-1:1), xc, yc
      integer ii,jj

      logical exceeds_th, transport_exceeds_th

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1)
      qmax = q(mq,1,1)
      do j = 1,my
         do i = 1,mx
             xc = xlower + (i-0.5)*dx
             yc = ylower + (j-0.5)*dy
             qmin = min(qmin,q(mq,i,j))
             qmax = max(qmax,q(mq,i,j))
              do ii = -1,1
                  do jj = -1,1
                        quad(ii,jj) = q(mq,i+ii,j+jj)
                  end do
              end do
              exceeds_th = transport_exceeds_th(blockno,
     &                     q(mq,i,j),qmin,qmax,quad, dx,dy,xc,yc, 
     &                     tag_threshold)

              if (exceeds_th) then
                  tag_patch = 1
                  return
              endif
         enddo
      enddo

      end
