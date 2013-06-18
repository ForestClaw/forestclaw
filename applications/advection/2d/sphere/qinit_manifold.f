      subroutine qinit_manifold(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux,
     &      xp,yp,zp)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux,this_block_idx
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, m
      double precision x,y,z

      do m = 1,meqn
         do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
               x = xp(i,j)
               y = yp(i,j)
               z = zp(i,j)

               if (x .le. 0) then
                  q(i,j,m) = 1.d0
               else
                  q(i,j,m) = 0.d0
               endif
            enddo
         enddo
      enddo

      return
      end
