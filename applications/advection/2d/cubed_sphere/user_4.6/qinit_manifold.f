      subroutine qinit_manifold(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux, xp,yp,zp)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux,this_block_idx
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision x,y,z, xlow, ylow, w

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(xlow,ylow,dx,dy,w)
            q(i,j,1) = w
         enddo
      enddo

      return
      end


      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer*8 cont, get_context
      integer blockno, get_block

      cont = get_context()
      blockno = get_block()

c      call mapc2m(xc,yc,xp,yp,zp)

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)


      fdisc = -xp
      end
