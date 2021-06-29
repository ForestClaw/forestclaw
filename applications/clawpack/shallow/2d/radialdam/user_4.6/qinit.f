      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Shallow water with radial dam break problem, h = hin inside
c     # circle specified in fdisc.f
c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision xlow, ylow,win
      integer i,j

      integer blockno, fc2d_clawpack46_get_block

      double precision hin, hout
      common /comic/ hin,hout

      blockno = fc2d_clawpack46_get_block()

      do i=1-mbc,mx+mbc
         xlow = xlower + (i-1.d0)*dx
         do j=1-mbc,my+mbc
            ylow = ylower + (j-1.d0)*dy
            call cellave2(blockno,xlow,ylow,dx,dy,win)
            q(i,j,1) = hin*win + hout*(1.d0-win)
            q(i,j,2) = 0.d0
            q(i,j,3) = 0.d0
            if (meqn .eq. 4) then
               q(i,j,4) = 0.d0
            endif
         end do
      end do
      return
      end
