      subroutine tag4refinement(mx,my,mbc,meqn,
     &      xlower,ylower, dx,dy,q,initflag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, initflag, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qmin, qmax, qavg

      end

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,qcoarsened, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, maux
      double precision xlower, ylower, dx, dy
      double precision qcoarsened(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision xc,yc,bmount, s

      integer i,j, mq

      double precision ampl, tol1, tol2, eta

      common /cprob/ ampl, tol1, tol2


c     # The difference between this and the true "refinement" above is
c     # that we can't check ghost cells here.  Also, we may make the
c     # coarsening criteria different from the refinement criteria.
c     # Also, we don't check for an init_flag, since it is unlikely that
c     # we would coarsen an initial grid.

      tag_patch = 0
      do i = 1,mx
         do j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (i-0.5)*dy
            s = bmount(xc,yc)
            eta = qcoarsened(i,j,1) + s
            if (abs(eta) .gt. tol1) then
               tag_patch = 1
               return
            endif
         enddo
      enddo



      end
