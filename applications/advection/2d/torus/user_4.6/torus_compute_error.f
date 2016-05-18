      subroutine torus_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t, q,error)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,m
      double precision xc,yc, qexact

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            error(i,j,1) = abs(q(i,j,1) - qexact(blockno,xc,yc,t));
         enddo
      enddo


      end


      double precision function qexact(blockno,xc,yc,t)
      implicit none

      integer blockno
      double precision xc,yc,t
      double precision x0, y0, u0, v0
      double precision q0

      u0 = 0.d0
      v0 = 0.d0

c     # Assume velocity is horizontal;  unit speed.
      qexact = q0(blockno, xc-t,yc)


      end
