       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer i,j
      double precision xc,yc, xlow, ylow,w,xp,yp,zp
      double precision q0_physical
      integer blockno

      blockno = 0
      do j = 1-mbc,my+mbc
          do i = 1-mbc,mx+mbc     
              xc = xlower + (i-0.5)*dx
              yc = ylower + (j-0.5)*dy
              if (initchoice .eq. 0) then
c                 # Discontinuous solution
                  xlow = xlower + (i-1)*dx
                  ylow = ylower + (j-1)*dy
                  call cellave2(blockno,xlow,ylow,dx,dy,w)
                  q(1,i,j) = w
              elseif (initchoice .ge. 1) then               
                  call mapc2m_cylinder(xc,yc,xp,yp,zp)
                  q(1,i,j) = q0_physical(xp,yp,zp)
              endif
          enddo
      enddo

      return
      end
