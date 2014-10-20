      subroutine swirl_tag4refinement(mx,my,mbc,meqn,blockno,
     &      xlower,ylower,dx,dy,q,init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag, blockno
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qmin, qmax
      double precision dq, dqi, dqj
      double precision xp,yp,zp
      integer*8 cont, get_context

      cont = get_context()

      qmin = 100.d0
      qmax = -100.d0
      tag_patch = 0


c     # Refine based only on first variable in system.
      mq = 1
      do i = 1,mx
         do j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
c            if (xp .ge. 0.5d0) then
c               tag_patch = 1
c               return
c            endif

            if (init_flag .eq. 1 .and. .false.) then
               if (abs(xc-0.5) .lt. dy) then
                  tag_patch = 1
                  return
               endif
            else
               qmin = min(q(i,j,mq),qmin)
               qmax = max(q(i,j,mq),qmax)
               if (qmax - qmin .gt. 0.25d0) then
                  tag_patch = 1
                  return
               endif
            endif
         enddo
      enddo

      end

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine swirl_tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, qcoarsened, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      double precision xlower, ylower, dx, dy
      double precision qcoarsened(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

c     # The difference between this and the true "refinement" above is
c     # that we can't check ghost cells here.  Also, we may make the
c     # coarsening criteria different from the refinement criteria.
c     # Also, we don't check for an init_flag, since it is unlikely that
c     # we would coarsen an initial grid.

      tag_patch = 0
      qmin = 100.d0
      qmax = -100.d0
      mq = 1
      do i = 1,mx
         do j = 1,my
            qmin = min(qcoarsened(i,j,mq),qmin)
            qmax = max(qcoarsened(i,j,mq),qmax)
            if (qmax - qmin .gt. 0.5d0) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
