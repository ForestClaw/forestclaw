      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno


      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      integer example
      common /example_comm/ example  

      double precision init_radius
      common /initradius_comm/ init_radius

      integer*8 cont, get_context
      double precision th, r0,x0,y0,z0,r
      double precision xc1, yc1

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      r0 = init_radius
      if (example .eq. 0 .or. example .eq. 1) then
          th = pi2*(0.25 + 1.d0/32.d0)
          x0 = cos(th)
          y0 = sin(th)
          z0 = alpha
      elseif (example .eq. 2) then
          xc1 = 0.25
          yc1 = 0.125
          call mapc2m_torus(xc1,yc1,x0,y0,z0,alpha,beta)
      endif

c     # Distance from thc
      r = sqrt((xp-x0)**2 + (yp-y0)**2 + (zp-z0)**2)

      fdisc = r - r0 

      end
