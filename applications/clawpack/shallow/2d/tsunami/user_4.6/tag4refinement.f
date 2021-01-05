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
      double precision qmin, qmax, xc, yc, bathy, hij
      


      tag_patch = 0

      cont = get_context()

      bathy = -1

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq) + bathy
      qmax = q(1,1,mq) + bathy
      do j = 1,my
         do i = 1,mx
            hij = q(i,j,1)
            qmin = min(hij + bathy,qmin)
            qmax = max(hij + bathy,qmax)
            if ((qmax) .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
