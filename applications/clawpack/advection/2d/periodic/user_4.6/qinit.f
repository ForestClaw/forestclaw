      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none
      integer meqn, mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer initial_condition
      common /com_initcond/ initial_condition

      integer i, j, mq
      double precision xlow,ylow, wl
      integer blockno, fc2d_clawpack46_get_block

      blockno = fc2d_clawpack46_get_block()

      do mq = 1,meqn
          do i = 1-mbc,mx+mbc
              xlow = xlower + (i-1)*dx
              do j = 1-mbc,my+mbc
                  ylow = ylower + (j-1)*dy
                  if (initial_condition .eq. 0) then
                      call cellave2(blockno,xlow,ylow,dx,dy,wl)
                      q(i,j,1) = wl
                  else if (initial_condition .eq. 1) then
                       q(i,j,1) = exp(-30*((xlow+dx/2)**2 + 
     &                     (ylow+dy/2)**2))
                  endif
              enddo
          enddo
      enddo

      return
      end



