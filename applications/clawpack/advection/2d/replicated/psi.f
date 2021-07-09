      double precision function psi(xc,yc)
      implicit none

      double precision xc, yc

      double precision uvel,vvel, revs_per_s
      common /comm_velocity/ uvel,vvel,revs_per_s

      psi = -revs_per_s*(uvel*xc - vvel*yc)

      end
