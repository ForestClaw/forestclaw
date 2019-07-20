      subroutine mapc2m_annulus(ac,bc,xp,yp,zp)
      implicit none

      double precision ac,bc, xc,yc,xp,yp,zp, r
      double precision l1(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      double precision t1, t2, t

      xc = ac
      yc = bc

      t1 = theta(1)
      t2 = theta(2)

      r = beta + (1-beta)*yc
      t = 2*pi*(t1 + (t2-t1)*xc)
      xp = r*cos(t)
      yp = r*sin(t)
      zp = 0


      end
