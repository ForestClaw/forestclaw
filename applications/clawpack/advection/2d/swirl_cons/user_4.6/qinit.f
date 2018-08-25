      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i, j, mq, blockno, fc2d_clawpack46_get_block
      double precision xc,yc, xp, yp, zp, xlow, ylow, w, tol
      double precision dxc,xm,ym

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      double precision pi
      integer example

      common /compi/ pi
      common /comex/ example

      cont = get_context()

      blockno = fc2d_clawpack46_get_block()

      do mq = 1,meqn
          do i = 1-mbc,mx+mbc
              xc = xlower + (i-0.5d0)*dx
              xlow = xlower + (i-1)*dx
              do j = 1-mbc,my+mbc
                  yc = ylower + (j-0.5d0)*dy
                  ylow = ylower + (j-1)*dy
                  if (example .le. 2) then
                      q(i,j,1) = 1                      
                  elseif (example .eq. 3) then
                      call cellave2(blockno,xlow,ylow,dx,dy,w)
                      q(i,j,1) = w
                  elseif (example .eq. 4) then
                      if (fclaw2d_map_is_used(cont)) then
                          call fclaw2d_map_c2m(cont,
     &                           blockno,xc,yc,xp,yp,zp)
                      else
                          xp = xc
                          yp = yc
                      endif
                     dxc = 0.25d0
                     xm = 0.5-dxc/2
                     ym = 0.5
                     w = dxc/2.d0
                     q(i,j,mq) = 0
                     if (abs(xp-xm) .le. w) then
                         q(i,j,mq) = 1
                     endif
                  endif
              enddo
          enddo
      enddo

      return
      end
