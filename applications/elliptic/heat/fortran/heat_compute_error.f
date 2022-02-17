      subroutine heat_compute_error(blockno, mx,my,mbc,mfields,
     &      dx,dy,xlower,ylower,t,rhs,error,soln)
      implicit none

      integer mx,my,mbc,mfields, blockno
      double precision dx, dy, xlower, ylower, t
      double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

      integer*8 cont, get_context

      integer i,j,m
      double precision xc,yc, heat_qexact

      cont = get_context()

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
            yc = ylower + (j-0.5)*dy
            xc = xlower + (i-0.5)*dx

            soln(i,j,1) = heat_qexact(xc,yc)
            do m = 1,mfields
               error(i,j,m) = rhs(i,j,m) - soln(i,j,1)
            end do
         enddo
      enddo


      end
