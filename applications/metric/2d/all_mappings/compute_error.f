      subroutine compute_error(blockno,mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t,q,error)
     &      bind(c,name="compute_error")
      implicit none

      integer meqn,mbc,mx,my
      double precision xlower,ylower,dx,dy, t
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision xp,yp,zp, alpha,kappa
      double precision xc,yc
      integer*8 cont, get_context
      integer blockno
      logical fclaw2d_map_is_torus
      logical fclaw2d_map_is_sphere
      integer i,j

c     # Copied from 'fclaw2d_metric_terms.i'

      double precision beta
      common /comtorus/ beta

      double precision pi
      common /compi/ pi

      cont = get_context()

      do i = 1,my
         do j = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # Error in Jacobian, i.e. area(i,j)
            error(i,j,1) = 0  !! Not yet figured out

c           # Error in curvature
            if (fclaw2d_map_is_sphere(cont)) then
c              # --user:example=6 (pillowdisk); --user:example=7 (cubedsphere)
               error(i,j,2) = abs(q(i,j,2)-1.d0)
            elseif (fclaw2d_map_is_torus(cont)) then
c              # --user:example=8 (torus)
               kappa = (1 + 2*beta*cos(2*pi*yc))/
     &               (2*beta*(1 + beta*cos(2*pi*yc)))
               error(i,j,2) = abs(q(i,j,2) - kappa)
            else
c              # --user:example=0-5  (flat;  curvature = 0)
               error(i,j,2) = 0
            endif
         enddo
      enddo

      end
