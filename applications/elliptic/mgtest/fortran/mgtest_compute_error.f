      subroutine mgtest_compute_error(blockno, mx,my,mbc,mfields,
     &      dx,dy,xlower,ylower,t,rhs,error,soln)
      implicit none

      integer mx,my,mbc,mfields, blockno
      double precision dx, dy, xlower, ylower, t
      double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

      integer*8 cont, get_context

      integer i,j,m
      double precision xc,yc, mgtest_qexact
      double precision xc1, yc1, zc1

      cont = get_context()

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
            yc = ylower + (j-0.5)*dy
            xc = xlower + (i-0.5)*dx

            soln(i,j,1) = mgtest_qexact(xc,yc)
            error(i,j,1) = rhs(i,j,1) - soln(i,j,1)
         enddo
      enddo


      end
