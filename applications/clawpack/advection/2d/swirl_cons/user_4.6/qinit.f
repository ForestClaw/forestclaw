      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i, j, mq, blockno, fc2d_clawpack46_get_block
      double precision xc,yc, xp, yp, zp, xlow, ylow, w, tol
      double precision dxc,xm

      double precision pi
      integer example

      common /compi/ pi
      common /comex/ example

      integer*8 map_context_ptr, get_context

      map_context_ptr = get_context()

      blockno = fc2d_clawpack46_get_block()

      do mq = 1,meqn
          do i = 1-mbc,mx+mbc
              xc = xlower + (i-0.5d0)*dx
              xlow = xlower + (i-1)*dx
              do j = 1-mbc,my+mbc
                  yc = ylower + (j-0.5d0)*dy
                  ylow = ylower + (j-1)*dy
                  q(i,j,mq) = 1
                  if (example .eq. 3) then
                      dxc = 0.25d0
                      xm = 0.5-dxc/2
                      w = dxc/2.d0
                      if (abs(yc-xm) .le. w) then
                          q(i,j,mq) = (xm-yc)/(w)
c                           q(i,j,1) = 1.d0/(mx*0.125*0.125)
                      endif 
                  endif
              enddo
          enddo
      enddo

      return
      end
