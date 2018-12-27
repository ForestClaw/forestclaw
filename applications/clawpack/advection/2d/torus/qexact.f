c     !! This is used to get the error
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
      double precision Hsmooth

      cont = get_context()

      xc = xc1
      yc = yc1
      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)


c     # Sphere centered at (1,0,r0) on torus
      r0 = 0.4d0
      r = sqrt((xp - 1.0)**2 + yp**2 + (zp-r0)**2)
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      end

      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end
