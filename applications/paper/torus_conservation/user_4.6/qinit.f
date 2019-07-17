      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j

c     # This is not the torus mapping, but rather maps the brick to
c     # a unit square      


      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            q(i,j,1) = 1.d0
         enddo
      enddo

      return
      end
