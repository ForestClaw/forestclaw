      subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer meqn, mbc, mx, my, maux,maxmx,maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      double precision x,y,xp,yp,zp
      double precision gaussian_sum
      double precision r,hmax,b,c
      integer blockno, fc2d_clawpack46_get_block
      integer*8 cont, get_context

      cont = get_context()
      blockno = fc2d_clawpack46_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            x = xlower + (i-1)*dx
            y = ylower + (j-1)*dy
            call fclaw_map_2d_c2m(cont,blockno,x,y,xp,yp,zp)
            q(i,j,1) = gaussian_sum(xp,yp,zp) + 0.2d0
            q(i,j,2) = 0
            q(i,j,3) = 0
            q(i,j,4) = 0
         enddo
      enddo

      return
      end
