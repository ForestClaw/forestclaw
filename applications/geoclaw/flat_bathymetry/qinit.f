c
c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
      use geoclaw_module, only: sea_level, coordinate_system
      use geoclaw_module, only: grav

      implicit none
      integer i,j,meqn,mbc,mx,my,maux
      double precision xlower,ylower,dx,dy,xi,yj
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      double precision z, ze

c     # It would be nice to get a simple wave here, i.e. a single right 
c     # going wave

      q = 0.d0
      do i = 1-mbc,mx+mbc
c         # (xi, yj) is the location of the physical domain
          xi = xlower + (i-0.5d0)*dx
          do j=1-mbc,my+mbc
              yj = ylower + (j-0.5d0)*dy
              q(1,i,j) = sea_level - aux(1,i,j)
              q(2,i,j) = 0
              q(3,i,j) = 0
          enddo
      enddo

      return
      end
