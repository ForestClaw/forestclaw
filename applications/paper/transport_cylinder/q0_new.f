c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision r_cyl, h_cyl
      common /cylinder_comm/ r_cyl, h_cyl

      double precision xc0, yc0, r0
      common /cylinder_init_comm/ xc0, yc0, r0

      double precision xp0, yp0, zp0
      double precision q0, Hsmooth, r
      integer k

      call mapc2m_cylinder(xc0,yc0,xp0,yp0,zp0)

c     # Sphere centered at (0.5,0.5,0) on swirl
      if (initchoice .eq. 1) then
          q0 = 1.d0
      elseif (initchoice .eq. 2) then
          r = sqrt((xp - xp0)**2 + (yp-yp0)**2 + (zp-zp0)**2)
          q0 = 1 - Hsmooth(r - r0)
      endif
      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end

      double precision function q0_init(xc,yc)
      implicit none 

      double precision xc,yc

      double precision xp, yp, zp
      double precision q0_physical

      call mapc2m_cylinder(xc,yc,xp,yp,zp)
      q0_init = q0_physical(xp,yp,zp)

      end
      







