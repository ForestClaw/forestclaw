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

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision x0, y0, z0
      double precision q0, Hsmooth, r, r0, th
      integer k

c     # Sphere centered at (0.5,0.5,0) on swirl
      r0 = init_radius
      if (initchoice .eq. 1) then
          q0 = 1.d0
      elseif (initchoice .eq. 2) then
          x0 = 1
          y0 = 0
          z0 = alpha

          r = sqrt((xp - x0)**2 + (yp-y0)**2 + (zp-z0)**2)
          !!q0 = Hsmooth(r + r0) - Hsmooth(r - r0)
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

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision xp, yp, zp
      double precision q0_physical

      call mapc2m_torus2(xc,yc,xp,yp,zp)
      q0_init = q0_physical(xp,yp,zp)

      end
      







