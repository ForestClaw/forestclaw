      double precision function bmount(xc,yc)
c     ============================================================
c
      implicit double precision (a-h,o-z)
c
      common /comsphere/ Rsphere, Omega
      common /comic/ Px,Py,Pz


      pi = 4.d0*datan(1.d0)
      bmount = -4.d4

c     # compute latitude:
      call mapc2m(xc,yc,xp,yp,zp)
      theta = dasin((xp*Px + yp*Py + zp*Pz) / Rsphere)

      thetam = -pi/6.d0 !position of ridge
      bmount = bmount + 3.d4*dexp(-1000.d0*(theta-thetam)**2)
      return
      end
