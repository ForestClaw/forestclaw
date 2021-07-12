      subroutine mapc2m_annulus2(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      double precision theta, r

      call map_comp2annulus(xc,yc,theta,r)

      xp = r*cos(theta)
      yp = r*sin(theta)
      zp = 0

      end
