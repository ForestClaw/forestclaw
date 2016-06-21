      subroutine torus_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t, q,error)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,m
      double precision xc,yc, qexact

      integer ex_comm, example
      common /comm_example/ ex_comm

      example = ex_comm


c     # Assume a single field variable only
      do j = 1,my
         yc = ylower + (j-0.5)*dy
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            error(i,j,1) = q(i,j,1) - qexact(blockno,xc,yc,t);
         enddo
      enddo


      end


      double precision function qexact(blockno,xc,yc,t)
      implicit none

      integer blockno
      double precision xc,yc,t
      double precision x0, y0, u0, v0
      double precision q0,qc
      double precision xc1, yc1, th,xp,yp,zp

      integer*8 cont, get_context

      double precision u0_comm,v0_comm,revs_comm
      common /comm_velocity/ u0_comm,v0_comm, revs_comm

      integer ex_comm, example
      common /comm_example/ ex_comm

      double precision pi
      common /compi/ pi

      double precision revs_per_sec_comm, r0_comm
      common /comm_torus/ revs_per_sec_comm, r0_comm

      example = ex_comm

      if (example .eq. 6) then
         u0 = revs_comm*u0_comm
         v0 = revs_comm*v0_comm

c        # Assume velocity is horizontal;  unit speed.
         qc = q0(blockno, xc - u0*t,yc - v0*t, t)
      elseif (example .eq. 7) then
c        # rigid body rotation on a torus
         qc = q0(blockno,xc,yc,t)

      endif
      qexact = qc


      end
