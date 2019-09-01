      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp, r2, w1, w2, thc, g, r

      double precision pi
      common /compi/ pi

      double precision alpha, revs_per_s
      common /torus_comm/ alpha, revs_per_s


      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      thc = pi/2.0
      thc = 0
      w1 = pi/8.d0
      w2 = pi/4.0

c     # atan2(y,x) \in [-pi,pi]      
      th = atan2(yp,xp)

c     # asin(zp/alpha) \in [-pi/2,pi/2]      
      g = asin(zp/alpha)

c     # Distance from thc
      r2 = (th-thc)**2

      fdisc = r2 - w1**2 

      end
