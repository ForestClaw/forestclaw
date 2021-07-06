c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision x0(2), y0(2), z0(2)
      common /qinit_comm/ x0, y0, z0

      integer example
      common /example_comm/ example

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision b_init, c_init
      common /init_parms/ b_init, c_init


      double precision q0, r0, r, cb(2)
      double precision theta, th, thc(2)
      double precision phi, l, lc(2)
      double precision c,b,cbell, hi, hmax, rarg
      double precision Hsmooth
      integer k

c     # Sphere centered at (0.5,0.5,0) on swirl
      if (initchoice .eq. 1) then
          q0 = 1.d0
      elseif (initchoice .eq. 2) then
          b = b_init
          c = c_init
          x0(1) = cos(pi/6.d0)
          x0(2) = cos(pi/6.d0)
          y0(1) = 0.5d0
          y0(2) = -0.5d0
          z0(1) = 0
          z0(2) = 0

          r0 = 0.35
          q0 = 0
          do  k = 1,2
              r = sqrt((xp - x0(k))**2 + (yp-y0(k))**2 
     &                 + (zp-z0(k))**2)
              q0 = q0 + 1 - Hsmooth(r - r0)
          end do
          q0 = b + c*q0
      elseif (initchoice .ge. 3) then
c         # Cosine Bells        
          hmax = 1
          b = b_init
          c = c_init
          r0 = 0.5d0

          call map2spherical(xp,yp,zp,theta,phi)
          l = theta
          th = phi
          
          if (example .eq. 1) then
              lc(1) = pi
              lc(2) = pi
              thc(1) = pi/3.d0
              thc(2) = -pi/3.d0
          elseif (example .eq. 2 .or. example .eq. 4) then
              lc(1) = 5*pi/6.d0
              lc(2) = 7*pi/6.d0
              thc(1) = 0
              thc(2) = 0
          elseif (example .eq. 3) then
              lc(1) = 3*pi/4.d0
              lc(2) = 5*pi/4.d0
              thc(1) = 0
              thc(2) = 0
          endif

          cbell = b
          do k = 1,2
              rarg = sin(thc(k))*sin(th) + 
     &                    cos(thc(k))*cos(th)*cos(l-lc(k))
              if (abs(rarg) > 1) then
                   write(6,*)  'rarg > 1'
                   stop
              endif

              r = acos(rarg)
              hi = hmax/2*(1 + cos(pi*r/r0))
              if (r .le. r0) then
                  cbell = b + c*hi
              endif
          end do
          q0 = cbell

      endif
      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none      

      double precision r

      double precision sharpness
      common /hsmooth_parms/ sharpness

      Hsmooth = (tanh(r/sharpness) + 1)/2.d0

      end

      double precision function q0_init(xc,yc)
      implicit none 

      double precision xc,yc

      integer example
      common /example_comm/ example

      double precision xp, yp, zp, q0
      double precision q0_physical

      call mapc2m_spherical(xc,yc,xp,yp,zp)

      q0_init = q0_physical(xp,yp,zp)

      end
      







