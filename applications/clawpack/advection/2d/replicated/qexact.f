      double precision function qexact(blockno,xc,yc,t)
      implicit none

      integer blockno
      double precision xc,yc,t
      double precision x0, y0, u0, v0
      double precision q0,qc

      double precision u0_comm,v0_comm,revs_comm
      common /comm_velocity/ u0_comm,v0_comm, revs_comm

      u0 = revs_comm*u0_comm
      v0 = revs_comm*v0_comm

c     # Assume velocity is horizontal;  unit speed.
      qc = q0(blockno, xc - u0*t,yc - v0*t)

      qexact = qc

      end



      double precision function  q0(blockno,xc1,yc1)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      double precision xc1, yc1
      integer blockno
      integer*8 cont, get_context
      double precision r,r0
      logical fclaw_map_is_used
      double precision Hsmooth,h1,h2

      double precision xloc(0:4),yloc(0:4),zloc(0:4)

      double precision pi
      common /compi/ pi

      data xloc /0, 1, 1, 0, 0.5d0/
      data yloc /0, 0, 1, 1, 0.5d0/

      cont = get_context()

      if (xc1 .lt. 0) then
         xp = 1 - (-xc1 - int(-xc1))
      else
         xp = xc1 - int(xc1)
      endif
      if (yc1 .lt. 0) then
         yp = 1-(-yc1 - int(-yc1))
      else
         yp = yc1 - int(yc1)
      endif
      r0 = 0.2d0
      rp = sqrt((xp-xloc(4))**2 + (yp-yloc(4))**2)
      h1 = Hsmooth(rp+r0)
      h2 = Hsmooth(rp-r0)
      q0 = h1 - h2

      end

      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end
