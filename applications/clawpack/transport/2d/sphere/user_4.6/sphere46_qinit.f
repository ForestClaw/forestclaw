      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i, j, blockno, fc2d_clawpack46_get_block
      double precision xc,yc, xp, yp, zp, xlow, ylow, w

      double precision q0_physical

      integer*8 cont, fclaw_map_get_context

      integer initchoice
      common /initchoice_comm/ initchoice

      cont = fclaw_map_get_context()

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
                  q(i,j,1) = 1.d0
              else
                  call fclaw_map_2d_c2m(cont,blockno,xc,yc,xp,yp,zp)                  
                  q(i,j,1) = q0_physical(xp,yp,zp)                  
              endif
          enddo
      enddo

      return
      end


      
