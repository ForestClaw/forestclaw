      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision xc,yc,xp,yp,zp
      double precision gaussian_sum

      integer*8 cont, fclaw_map_get_context
      integer blockno, fc2d_clawpack46_get_block

      cont = fclaw_map_get_context()
      blockno = fc2d_clawpack46_get_block()

      do j = 1-mbc,my+mbc
         yc = ylower + (j-0.5d0)*dy
         do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            call fclaw_map_2d_c2m(cont,
     &            blockno,xc,yc,xp,yp,zp)
            q(i,j,1) = gaussian_sum(xp,yp,zp)
         enddo
      enddo

      return
      end
