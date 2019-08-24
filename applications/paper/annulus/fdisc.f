      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example      

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision x0, y0, r0, r, rc,th

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      if (example .ne. 2) then
c         # Annular or horizontal flow            
          rc = (1 + beta)/2.d0
          th = pi2*(0.25 + 1.d0/32.d0)
          x0 = rc*cos(th)
          y0 = rc*sin(th)
      elseif (example .eq. 2) then
c         # Vertical flow              
          th = pi/2
          rc = beta + (1-beta)*0.625d0
          x0 = rc*cos(th)
          y0 = rc*sin(th)
      endif

      r = sqrt((xp - x0)**2 + (yp-y0)**2)
      r0 = init_radius
      fdisc = r - r0

      end
