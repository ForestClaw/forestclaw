      subroutine tag4refinement(mx,my,mbc,
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

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      


      tag_patch = 0

      cont = get_context()

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
c            if (abs(q(i,j,mq)) .gt. tag_threshold) then
            if (xc > 0.5) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
