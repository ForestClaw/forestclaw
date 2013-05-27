      subroutine qinit_manifold(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux,
     &      xp,yp,zp,xd,yd,zd,this_block_idx)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux,this_block_idx
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      integer i,j
      double precision x,y,z

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            x = xp(i,j)
            y = yp(i,j)
            z = zp(i,j)

            if (y .le. 0) then
               q(i,j,1) = 1.d0
            else
               q(i,j,1) = 0.d0
            endif
         enddo
      enddo

      return
      end
