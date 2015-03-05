      subroutine tag4refinement(mx,my,mbc,meqn,  xlower,ylower,
     &      dx,dy,blockno, q, refine_threshold,init_flag, tag_patch)
     &
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer init_flag, blockno
      double precision refine_height, refine_threshold
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc
      double precision xp,yp,zp

      double precision qmin,qmax
      double precision zm, mountain_height
      integer*8 cont, get_context

      cont = get_context()

      tag_patch = 0
      refine_height = 200

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,1)
      qmax = q(1,1,1)
      do j = 1,my
         do i = 1,mx
            qmin = min(qmin,q(i,j,1))
            qmax = max(qmax,q(i,j,1))
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
            zm = mountain_height(xp)
            if ((zm .lt. yp .and. yp .lt. zm + refine_height) .or.
     &            qmax - qmin .gt. refine_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
