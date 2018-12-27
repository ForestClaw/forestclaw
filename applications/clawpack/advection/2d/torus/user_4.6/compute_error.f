      subroutine torus46_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t, q,error)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,m
      double precision xc,yc, qexact

      integer ex_comm, example
      common /comm_example/ ex_comm

      example = ex_comm


c     # Assume a single field variable only
      do j = 1,my
         yc = ylower + (j-0.5)*dy
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            error(i,j,1) = q(i,j,1) - qexact(blockno,xc,yc,t);
         enddo
      enddo


      end
