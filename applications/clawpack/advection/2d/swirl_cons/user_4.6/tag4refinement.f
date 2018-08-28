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
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (q(i,j,mq) .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
            
c            xc = xlower + (i-0.5)*dx
c            yc = ylower + (j-0.5)*dy
c            if (fclaw2d_map_is_used(cont)) then
c               call fclaw2d_map_c2m(cont,
c     &         blockno,xc,yc,xp,yp,zp)
c            else
c               xp = xc
c               yp = yc
c            endif
c
c            if (mod(blockno,2) .eq. 1) then
c               tag_patch = 1
c               return
c            endif

         enddo
      enddo

      end
