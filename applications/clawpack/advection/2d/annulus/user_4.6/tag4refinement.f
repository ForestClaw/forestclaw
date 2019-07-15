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

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer mapping
      common /mapping_comm/ mapping

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer i,j, mq
      double precision qmin, qmax, xc, yc, xc1, yc1, zc1
      double precision r, rrefine, x,y, rw, ravg, th
      logical constant_theta, constant_r



      tag_patch = 0

      cont = get_context()

      rrefine = (1 + beta)/2.d0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)

      ravg = (1 + beta)/2.d0
      rw = (1-beta)/4.d0    
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            r = beta + (1-beta)*yc1
            th = pi2*(theta(1) + (theta(2)-theta(1))*xc1)
c           # constant_theta = cos(th) .lt. 0 .and. abs(sin(th)) .lt. 0.5
            !! constant_theta = t1 .le. th .and. th .le. t2
            constant_theta = th .gt. pi/2.d0
            constant_r = r > ravg                
            if (refine_pattern .eq. 0 .and. constant_theta) then 
                tag_patch = 1
                return
            elseif (refine_pattern .eq. 1 .and. constant_r) then
                tag_patch = 1
                return
            endif

c            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)
c            call annulus_transform_coordinates(xc1,yc1,x,y,mapping)
c            r = beta + (1-beta)*y
c            qmin = min(q(i,j,mq),qmin)
c            qmax = max(q(i,j,mq),qmax)
c            if (q(i,j,mq) .gt. tag_threshold .and. 
c     &          q(i,j,mq) .lt. 1-tag_threshold) then
c                tag_patch = 1
c                return
cc                if (r .gt. rrefine) then
cc                   tag_patch = 1
cc                   return
cc                endif
c            endif
         enddo
      enddo
            
      end
