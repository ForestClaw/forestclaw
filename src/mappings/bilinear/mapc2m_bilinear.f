      subroutine mapc2m_bilinear(blockno,xc,yc,xp,yp,zp,quad)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno
      double precision quad(0:1,0:1,3)

      double precision a00(3), a01(3),a10(3), a11(3)
      integer m

c     # Maps [0,1]x[0,1] to a bilinear map with physical 
c     # corners at [quad(0,0,:), quad(0,1,:), quad(1,0,:), quad(1,1,:)]


      do m = 1,3
         a00(m) = quad(0,0,m)
         a01(m) = quad(1,0,m) - quad(0,0,m)
         a10(m) = quad(0,1,m) - quad(0,0,m)
         a11(m) = quad(1,1,m) - quad(1,0,m) -
     &         quad(0,1,m) + quad(0,0,m)
      enddo

c     #  T(xi,eta) = a00 + a01*xi + a10*eta + a11*xi*eta

      xp = a00(1) + a01(1)*xc + a10(1)*yc + a11(1)*xc*yc
      yp = a00(2) + a01(2)*xc + a10(2)*yc + a11(2)*xc*yc
      zp = a00(3) + a01(3)*xc + a10(3)*yc + a11(3)*xc*yc

      end
