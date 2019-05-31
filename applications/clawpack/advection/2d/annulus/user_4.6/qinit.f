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

c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision r,r0, x0, y0, z0, q0, ravg
      double precision Hsmooth

      double precision beta
      common /annulus_comm/ beta

      double precision init_radius
      common /initradius_comm/ init_radius

c     # Sphere centered at (x0,0,0) on annulus
c     # Outer radius  = 1; inner radius = beta
c     # average inner and outer radii to center sphere
      r0 = init_radius
      ravg = (1 + beta)/2.d0
      x0 = -ravg
      y0 = 0

      r = sqrt((xp - x0)**2 + (yp - y0)**2)
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end


