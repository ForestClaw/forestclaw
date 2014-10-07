      subroutine qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux,this_block_idx
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      integer blockno, clawpack_get_block
      double precision x,y,z, xlow, ylow, w

      blockno = clawpack_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,w)
            q(i,j,1) = w
         enddo
      enddo

      return
      end


      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)


      fdisc = -xp
      end
