      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      integer*8 cont, get_context
      integer blockno, fc2d_clawpack46_get_block
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision xc,yc,xp,yp,zp
      double precision cosine_bell_sum

      double precision a,b

      cont = get_context()
      blockno = fc2d_clawpack46_get_block()

      a = -0.8d0
      b = 0.9d0

      do j = 1-mbc,my+mbc
         yc = ylower + (j-0.5d0)*dy
         do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            call fclaw2d_map_c2m(cont,
     &            blockno,xc,yc,xp,yp,zp)
            q(i,j,1) = cosine_bell_sum(xp,yp,zp)
            if (meqn .eq. 2) then
c              # Set non-linear relationship between two tracers
               q(i,j,2) = a*q(i,j,1)**2 + b
            endif
         enddo
      enddo

      return
      end
