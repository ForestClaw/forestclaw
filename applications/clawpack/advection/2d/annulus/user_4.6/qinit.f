      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      integer blockno, fc2d_clawpack46_get_block
      double precision xlow, ylow, w, xc, yc
      double precision xp, yp, zp

      integer*8 cont, get_context

      double precision q0_physical

      integer initchoice
      common /initchoice_comm/ initchoice

      cont = get_context()

      blockno = fc2d_clawpack46_get_block()

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            if (initchoice .eq. 0) then
               xlow = xlower + (i-1)*dx
               ylow = ylower + (j-1)*dy
               call cellave2(blockno,xlow,ylow,dx,dy,w)
               q(i,j,1) = w
            elseif (initchoice .eq. 1) then
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)                  
               q(i,j,1) = q0_physical(xp,yp,zp)
            else
               q(i,j,1) = 1.d0
            endif
         enddo
      enddo

      return
      end
