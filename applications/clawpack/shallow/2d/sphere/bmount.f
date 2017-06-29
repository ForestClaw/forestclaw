      double precision function bmount(xc,yc)
      implicit none

      double precision xc,yc
      double precision Rsphere, Omega, Px, Py, Pz
      double precision theta,thetam
      double precision pi, xp, yp, zp

      common /comsphere/ Rsphere, Omega
      common /comic/ Px,Py,Pz
      common /compi/ pi

      bmount = -4.d4

c     # compute latitude:
      call mapc2m(xc,yc,xp,yp,zp)
      theta = asin((xp*Px + yp*Py + zp*Pz) / Rsphere)

      thetam = -pi/6.d0 !position of ridge
      bmount = bmount + 3.d4*exp(-1000.d0*(theta-thetam)**2)
      return
      end
