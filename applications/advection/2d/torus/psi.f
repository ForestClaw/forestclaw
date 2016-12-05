      double precision function psi(blockno,xc,yc,t)
      implicit none

      double precision xc, yc, t
      integer blockno
      double precision pi, alpha
      double precision revs_per_s
      integer*8 cont, get_context

      double precision xc1, yc1, zc1, pi2
      integer example

      common /compi/ pi
      common /toruscomm/ alpha
      common /excomm_example/ example

      cont = get_context()

c     # (xc,yc) in each block is in [0,1]x[0,1].  This call maps
c     # maps that point to subregion of this unit domain, according
c     # to what brick the point occupies.  For example, in a 2x2 brick
c     # arrangement, the point (0.5,0.5) in brick 0 gets mapped to (0.25
c     #(0.25,0.25).  The same point in brick 2 gets mapped to
c     # (0.25,0.75).

      call fclaw2d_map_brick2c(cont,
     &      blockno,xc,yc,xc1,yc1,zc1)

      revs_per_s = 0.5d0

      pi2 = 2*pi

      if (example .eq. 0) then
c        # Rigid body rotation
         psi = (pi2*revs_per_s)*alpha*
     &         (pi2*yc1 + alpha*sin(pi2*yc1))
      elseif (example .eq. 1) then
c        # Twisted torus stream function (to be used with usual torus map)
         psi = (pi2*revs_per_s)*alpha*
     &         (pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
      endif

c      psi = psi

      end



      subroutine get_vel_psi_comp(blockno,xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(2),xd2(2), ds, vn, psi,t
      integer blockno

      vn = (psi(blockno,xd1(1),xd1(2),t) -
     &      psi(blockno,xd2(1),xd2(2),t))/ds

      end
