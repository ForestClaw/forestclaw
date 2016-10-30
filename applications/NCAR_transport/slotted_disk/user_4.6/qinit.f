      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      integer blockno, fc2d_clawpack46_get_block
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision xlow,ylow,w
      double precision r,hmax,b,c

      blockno = fc2d_clawpack46_get_block()

      CALL get_td_sdisk_parms(r,hmax,b,c)

      do j = 1-mbc,my+mbc
         ylow = ylower + (j-1)*dy
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(i,j,1) = b + c*w
         enddo
      enddo

      return
      end
