      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsening_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy

c # Set as 'coarsening_threshold' in fclaw_options
      double precision refine_height, coarsening_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq, igrid
      double precision qmin, qmax, xc,yc,zm,xp,yp,zp
      double precision mountain_height

      integer*8 cont, get_context

      cont = get_context()
      refine_height = 200

c     # If any grid is near the terrain, never coarsen.
      tag_patch = 0
      do igrid = 0,3
         do i = 1,mx
            do j = 1,my
               xc = xlower(igrid) + (i-0.5)*dx
               yc = ylower(igrid) + (j-0.5)*dy
               call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
               zm = mountain_height(xp)
               if (zm .lt. yp .and. yp .lt. zm + refine_height) then
c                 # Don't coarsen
                  tag_patch = 0
                  return
               endif
            enddo
         enddo
      enddo

c     # We are not near the terrain surface, so we will coarsen if
c     # we are at an advection front
      qmin = q0(1,1,1)
      qmax = q0(1,1,1)
      call get_minmax(mx,my,mbc,meqn,q0,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q1,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q2,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q3,qmin,qmax)
      if (qmax - qmin .lt. coarsening_threshold) then
         tag_patch = 1
         return
      endif


      end

      subroutine get_minmax(mx,my,mbc,meqn,q,qmin,qmax)

      implicit none
      integer mx,my,mbc,meqn
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j,mq

      mq = 1
      do i = 1,mx
         do j = 1,my
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
         enddo
      enddo

      end
