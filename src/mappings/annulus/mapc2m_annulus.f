      subroutine mapc2m_annulus(blockno,xc,yc,xp,yp,zp,beta,theta)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp, beta, theta(2)

      double precision pi,pi2
      common /compi/ pi,pi2

      double precision t1, t2, t, r

      t1 = theta(1)
      t2 = theta(2)

      t = t1 + (t2-t1)*xc
      r = beta + (1-beta)*yc
      xp = r*cos(2*pi*t)
      yp = r*sin(2*pi*t)
      zp = 0

      end
