      subroutine clawpack46_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision xp,yp,zp

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      integer i,j, mq, ii, jj
      double precision qmin, qmax, xc, yc, quad(-1:1,-1:1)

      logical exceeds_th, radialdam_exceeds_th
      
      tag_patch = 0

      cont = get_context()

      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            do ii = -1,1
               do jj = -1,1
                  quad(ii,jj) = q(i+ii,j+jj,mq)
               end do
            end do
            exceeds_th = radialdam_exceeds_th(blockno,
     &             q(i,j,mq),qmin,qmax,quad, dx,dy,xc,yc,
     &             tag_threshold)
            if (exceeds_th) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
