      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      integer i,j
      double precision xlow, ylow, w, xc,yc, q0, t0

      integer blockno, fc2d_clawpack46_get_block

      integer initchoice
      common /initchoice_comm/ initchoice

      blockno = fc2d_clawpack46_get_block()

      t0 = 0.d0
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            if (initchoice .eq. 0) then
c              # Discontinuous solution
               xlow = xlower + (i-1)*dx
               ylow = ylower + (j-1)*dy
               call cellave2(blockno,xlow,ylow,dx,dy,w)
               q(i,j,1) = w
            elseif (initchoice .eq. 1) then
c              # Smooth solution for computing the error
c              # Assumes (xc,yc) has not been mapped to brick layout
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               q(i,j,1) = q0(blockno,xc,yc)
            endif
         enddo
      enddo

      return
      end
