      subroutine b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision tperiod, pi2, vt, xll,yll, psi, pi
      integer blockno, fc2d_clawpack46_get_block

      common /compi/ pi
c

      blockno = fc2d_clawpack46_get_block()

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(i,j,2) = -(psi(blockno,xll, yll+dy,time) -
     &            psi(blockno,xll,yll,time)) / dy
            aux(i,j,3) =  (psi(blockno,xll+dx, yll,time) -
     &            psi(blockno,xll,yll,time)) / dx

         enddo
      enddo

      return
      end
