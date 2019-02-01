      subroutine mapc2m_torus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision alpha, r


      double precision pi, pi2
      common /compi/ pi, pi2

      r = 1 + alpha*cos(pi2*yc)

      xp = r*cos(pi2*xc)
      yp = r*sin(pi2*xc)
      zp = alpha*sin(pi2*yc)

      end

      subroutine mapc2m_torus_invert(xp,yp,zp,xc1,yc1,alpha)
      implicit none

      double precision xp,yp,zp,xc1,yc1
      double precision alpha

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision r

      xc1 = atan2(yp,xp)
      if (xc1 .lt. 0) then
          xc1 = xc1 + pi2
      endif
      xc1 = xc1/pi2

      yc1 = asin(zp/alpha)

      r = sqrt(xp**2 + yp**2);
      if (r .gt. 1 .and. yc1 .lt. 0) then
c         # Quad IV
          yc1 = yc1 + pi2
      elseif (r .le. 1) then
c         # Quads II and III      
          yc1 = pi-yc1
      endif
      yc1 = yc1/pi2

      end
