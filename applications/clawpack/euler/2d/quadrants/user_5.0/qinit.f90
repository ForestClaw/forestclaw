SUBROUTINE clawpack5_qinit(meqn,mbc,mx,my, &
     xlower,ylower,dx,dy,q,maux,aux)
  !!     =====================================================
  !!
  !!     # Set initial conditions for q.
  !!
  !!      # Data is piecewise constant with 4 values in 4 quadrants
  !!      # 2D Riemann problem from Figure 4 of
  !!        @article{csr-col-glaz,
  !!          author="C. W. Schulz-Rinne and J. P. Collins and H. M. Glaz",
  !!          title="Numerical Solution of the {R}iemann Problem for
  !!                 Two-Dimensional Gas Dynamics",
  !!          journal="SIAM J. Sci. Comput.",
  !!          volume="14",
  !!          year="1993",
  !!          pages="1394-1414" }


  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DIMENSION q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DIMENSION rpp(4),rpr(4),rpu(4),rpv(4)
  COMMON /cparam/  gamma,gamma1


  !! # First quadrant:
  rpp(1) = 1.5d0
  rpr(1) = 1.5d0
  rpu(1) = 0.d0
  rpv(1) = 0.d0

  !! # Second quadrant:
  rpp(2) = 0.3d0
  rpr(2) = 0.532258064516129d0
  rpu(2) = 1.206045378311055d0
  rpv(2) = 0.0d0

  !! # Third quadrant:
  rpp(3) = 0.029032258064516d0
  rpr(3) = 0.137992831541219d0
  rpu(3) = 1.206045378311055d0
  rpv(3) = 1.206045378311055d0

  !! # Fourth quadrant:
  rpp(4) = 0.3d0
  rpr(4) = 0.532258064516129d0
  rpu(4) = 0.0d0
  rpv(4) = 1.206045378311055d0

!! # location of four corners:
  xs = .8d0
  ys = .8d0

  DO i = 1-mbc,mx+mbc
     xcell = xlower + (i-0.5d0)*dx
     DO j = 1-mbc,my+mbc
        ycell = ylower + (j-0.5d0)*dy
        IF (xcell.GE.xs .AND. ycell.GE.ys) iq = 1
        IF (xcell.LT.xs .AND. ycell.GE.ys) iq = 2
        IF (xcell.LT.xs .AND. ycell.LT.ys) iq = 3
        IF (xcell.GE.xs .AND. ycell.LT.ys) iq = 4
        q(1,i,j) = rpr(iq)
        q(2,i,j) = rpr(iq)*rpu(iq)
        q(3,i,j) = rpr(iq)*rpv(iq)
        q(4,i,j) = rpp(iq)/gamma1 + 0.5d0*rpr(iq)*(rpu(iq)**2 + &
             rpv(iq)**2)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack5_qinit
