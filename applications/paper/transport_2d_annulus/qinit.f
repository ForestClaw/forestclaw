       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer i,j, blockno
      double precision xlow, ylow, w

      blockno = 0

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(1,i,j) = w
         enddo
      enddo

      return
      end


