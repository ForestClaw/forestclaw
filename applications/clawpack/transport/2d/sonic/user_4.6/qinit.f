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

      double precision q0_physical

      integer*8 cont, get_context
      logical fclaw_map_is_used

      integer ii,jj

      integer initchoice
      common /initchoice_comm/ initchoice

      cont = get_context()

      blockno = fc2d_clawpack46_get_block()

      do i = 1-mbc,mx+mbc
          xc = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
              yc = ylower + (j-0.5d0)*dy

              if (initchoice .eq. 0) then
c                 # Discontinuous solution
                  xlow = xlower + (i-1)*dx
                  ylow = ylower + (j-1)*dy
                  call cellave2(blockno,xlow,ylow,dx,dy,w)
                  q(i,j,1) = w
              elseif (initchoice .eq. 1) then
                  call fclaw_map_2d_c2m(cont,blockno,xc,yc,xp,yp,zp)                  
                  q(i,j,1) = q0_physical(xp,yp,zp)
              elseif (initchoice .eq. 2) then
                  q(i,j,1) = 1.d0
              elseif (initchoice .eq. 3) then
                  if (fclaw_map_is_used(cont)) then
                      call fclaw_map_2d_c2m(cont,
     &                           blockno,xc,yc,xp,yp,zp)
                  else
                      xp = xc
                      yp = yc
                  endif
                  dxc = 1.d0/8.d0
                  xm = 0.25
                  ym = 0.25
                  w = dxc/2.d0
                  q(i,j,mq) = 0                     
                  if (blockno .eq.2) then
                      do jj=7,8
                          do ii=1,2
                              q(ii,jj,1) = 1.d0
                          enddo
                      enddo
                  endif
                  if (blockno .eq. 1) then
                      do ii=5,7
c                         q(ii,my,1) = 1.d0
                      enddo
                  endif
              endif
          enddo
      enddo

      return
      end
