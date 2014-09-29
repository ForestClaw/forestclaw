      subroutine compute_error(meqn,mbc,mx,my,cont,blockno,
     &      xlower,ylower,dx,dy,curvature,error)
     &      bind(c,name="compute_error")
      implicit none

      integer meqn,mbc,mx,my
      double precision xlower,ylower,dx,dy
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision xp,yp,zp, alpha,kappa,pi
      double precision xc,yc
      integer*8 cont
      integer blockno
      logical fclaw2d_map_is_torus
      logical fclaw2d_map_is_sphere
      integer i,j

c     # Copied from 'metric_terms.i' (which, unfortunately, also includes xp,yp,zp,...
      double precision   curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      common /compi/ pi

      alpha = 0.4d0
      do i = 1,my
         do j = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)

            if (fclaw2d_map_is_sphere(cont)) then
               error(i,j) = abs(curvature(i,j)-1)
            elseif (fclaw2d_map_is_torus(cont)) then
               kappa = (1 + 2*alpha*cos(2*pi*yc))/
     &               (2*alpha*(1 + alpha*cos(2*pi*yc)))
               error(i,j) = abs(curvature(i,j) - kappa)
            else
               error(i,j) = 0
            endif
         enddo
      enddo

      end
