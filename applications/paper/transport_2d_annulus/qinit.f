c      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
c     &      xlower,ylower,dx,dy,q,maux,aux)
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer i,j
      double precision xlow, ylow, w, xc,yc, q0, t0
      double precision xc1, yc1, zc1, xp,yp,zp
      double precision q0_physical
      integer blockno

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision pi, pi2
      common /compi/ pi, pi2

      integer j1, j2

      blockno = 0

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            if (initchoice .eq. 0) then
               xlow = xlower + (i-1)*dx
               ylow = ylower + (j-1)*dy
               call cellave2(blockno,xlow,ylow,dx,dy,w)
               q(1,i,j) = w
            elseif (initchoice .eq. 1) then
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               call mapc2m_annulus(xc,yc,xp,yp,zp)
               q(1,i,j) = q0_physical(xp,yp,zp)
            else
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               if (yc .ge. 0.25 .and. yc .le. 0.75) then
                  q(1,i,j) = 1.d0
               else
                  q(1,i,j) = 0
               endif
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

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision th

c     # Sphere centered at (x0,0,0) on annulus
c     # Outer radius  = 1; inner radius = beta
c     # average inner and outer radii to center sphere
      ravg = (1 + beta)/2.d0
      th = pi2*(0.25 + 1.d0/32.d0)
      x0 = ravg*cos(th)
      y0 = ravg*sin(th)

      r = sqrt((xp - x0)**2 + (yp - y0)**2)

      r0 = init_radius
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end

