      subroutine clawpack5_qinit(meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer i, j, blockno, fc2d_clawpack5_get_block
      double precision xc,yc, xp, yp, zp, xlow, ylow, w

      double precision q0_physical

      integer*8 cont, get_context

      integer initchoice
      common /initchoice_comm/ initchoice

      cont = get_context()

      blockno = fc2d_clawpack5_get_block()

      do i = 1-mbc,mx+mbc
          xc = xlower + (i-0.5d0)*dx
          do j = 1-mbc,my+mbc
              yc = ylower + (j-0.5d0)*dy

              if (initchoice .eq. 0) then
c                 # Discontinuous solution
                  xlow = xlower + (i-1)*dx
                  ylow = ylower + (j-1)*dy
                  call cellave2(blockno,xlow,ylow,dx,dy,w)
                  q(1,i,j) = w
              elseif (initchoice .eq. 1) then
                  q(1,i,j) = 1.d0
              else
                  call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)                  
                  q(1,i,j) = q0_physical(xp,yp,zp)                  
              endif
          enddo
      enddo

      return
      end


      
