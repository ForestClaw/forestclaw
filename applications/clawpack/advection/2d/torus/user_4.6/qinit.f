      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      double precision xlow, ylow, w, xc,yc, t0
      double precision xc1, yc1, zc1
      double precision xp, yp, zp
      double precision q0, q0_physical

      integer blockno, fc2d_clawpack46_get_block
      integer*8 cont, get_context

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision pi
      common /compi/ pi

      blockno = fc2d_clawpack46_get_block()

      cont = get_context()

c     # This is not the torus mapping, but rather maps the brick to
c     # a unit square      


      t0 = 0.d0
      do j = 1-mbc,my+mbc
         yc = ylower + (j-0.5)*dy
         do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5)*dx
            if (initchoice .eq. 0) then
c              # Discontinuous solution
               xlow = xlower + (i-1)*dx
               ylow = ylower + (j-1)*dy
               call cellave2(blockno,xlow,ylow,dx,dy,w)
               q(i,j,1) = w
            elseif (initchoice .eq. 1) then

               call fclaw2d_map_c2m(cont,
     &                blockno,xc,yc,xp,yp,zp)

               q(i,j,1) = q0_physical(xp,yp,zp)
            elseif (initchoice .eq. 2) then
                  q(i,j,1) = 1.d0
            endif
         enddo
      enddo

      return
      end
