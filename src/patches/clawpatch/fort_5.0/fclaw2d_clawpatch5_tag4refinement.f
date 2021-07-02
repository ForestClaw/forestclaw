      subroutine fclaw2d_clawpatch5_fort_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, mq, ii, jj
      double precision qmin, qmax, xc,yc,quad(-1:1,-1:1)
      logical fclaw2d_clawpatch_minmax_exceeds_th, exceeds_th

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1)
      qmax = q(mq,1,1)
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(mq,i,j),qmin)
            qmax = max(q(mq,i,j),qmax)
            do ii = -1,1
               do jj = -1,1
                  quad(ii,jj) = q(mq,i+ii,j+jj)
               end do
            end do
            exceeds_th = fclaw2d_clawpatch_minmax_exceeds_th(
     &             blockno, q(mq,i,j),qmin,qmax,quad, dx,dy,xc,yc,
     &             tag_threshold)
            if (exceeds_th) then
c              # Refine this patch               
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
