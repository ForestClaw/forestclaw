      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
c     # Data is piecewise constant with 4 values in 4 quadrants
c     # 2D Riemann problem from Figure 4 of
c        @article{csr-col-glaz,
c          author="C. W. Schulz-Rinne and J. P. Collins and H. M. Glaz",
c          title="Numerical Solution of the {R}iemann Problem for
c                 Two-Dimensional Gas Dynamics",
c          journal="SIAM J. Sci. Comput.",
c          volume="14",
c          year="1993",
c          pages="1394-1414" }

c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision gamma, gamma1
      common /cparam/  gamma,gamma1

      double precision rpp(4),rpr(4),rpu(4),rpv(4)
      double precision xcell, ycell, xs, ys
      integer i,j, iq

c
c     # First quadrant:
      rpp(1) = 1.5d0
      rpr(1) = 1.5d0
      rpu(1) = 0.d0
      rpv(1) = 0.d0
c
c     # Second quadrant:
      rpp(2) = 0.3d0
      rpr(2) = 0.532258064516129d0
      rpu(2) = 1.206045378311055d0
      rpv(2) = 0.0d0
c
c     # Third quadrant:
      rpp(3) = 0.029032258064516d0
      rpr(3) = 0.137992831541219d0
      rpu(3) = 1.206045378311055d0
      rpv(3) = 1.206045378311055d0
c
c     # Fourth quadrant:
      rpp(4) = 0.3d0
      rpr(4) = 0.532258064516129d0
      rpu(4) = 0.0d0
      rpv(4) = 1.206045378311055d0
c
c     # location of four corners:
      xs = .8d0
      ys = .8d0
c
      do i = 1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         do j = 1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
            if (xcell.ge.xs .and. ycell.ge.ys) iq = 1
            if (xcell.lt.xs .and. ycell.ge.ys) iq = 2
            if (xcell.lt.xs .and. ycell.lt.ys) iq = 3
            if (xcell.ge.xs .and. ycell.lt.ys) iq = 4
            q(i,j,1) = rpr(iq)
            q(i,j,2) = rpr(iq)*rpu(iq)
            q(i,j,3) = rpr(iq)*rpv(iq)
            q(i,j,4) = rpp(iq)/gamma1 + 0.5d0*rpr(iq)*(rpu(iq)**2 +
     &                                 rpv(iq)**2)
         end do
      end do
      return
      end
