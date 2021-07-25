SUBROUTINE sgn(meqn,mbc,mx,xlower,dx,q,maux,aux,sgn_soln)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, maux
  DOUBLE PRECISION xlower, dx, t, dt

  DOUBLE PRECISION  q(meqn,1-mbc:mx+mbc)
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc)

  DOUBLE PRECISION  sgn_soln(1-mbc:mx+mbc)

  DOUBLE PRECISION :: alpha, breaking
  COMMON /common_dispersive/ alpha, breaking

  DOUBLE PRECISION :: grav, dry_tolerance, sea_level
  COMMON /common_swe/ grav, dry_tolerance, sea_level

  !! Triangular matrix (dummy matrices)   
  DOUBLE PRECISION dl(mx), d(mx), du(mx), b(mx)
  DOUBLE PRECISION hl, hr, hc, hl3, hr3, dx2
  INTEGER nrhs, ldb, info

  INTEGER i, ibc, m
  DOUBLE PRECISION dxeta, dxh, dxzb, d2xzb, T10, TNNp1, a3, dxsq
  double precision dm1, d0, dp1, bterm, xc

  logical zero_out_linsys


  !! Create Tridiagonal matrix

  !! 1d version of SOLVER
  !!
  !!   -alpha/3*d/dx(h^3 d/dx D) + h*(alpha(d/dx eta d/dx b + h/2*d^2/dx^2 b) + 1)*D = b
  !!
  !! Discrete version : 
  !!
  !! (from Basilisk)
  !!     (-alpha/3.*(hr3*D.x[1,0] + hl3*D.x[-1,0] - (hr3 + hl3)*D.x[])/sq(Delta) +
  !!           hc*(alpha*(dxeta*dxzb + hc/2.*d2x(zb)) + 1.)*D.x[]
  !!
  !! Discretization
  !!
  !!     -alpha/3*(hl3*D_{i-1} - 2*(hr3+hl3)/2*D_{i} + hr3*D_{i+1})/dx^2 +
  !!          hc*(alpha*(dxeta*dxdzb + hc/2*dx2zb) + 1)*D_{i} = b
  !! 
  !! If bathy = -1, terms invovling dxdb are all 0, and we just get : 
  !!
  !!  -alpha/3*(hl3*D_{i-1}/dxsq +(hc - 2*(hr3+hl3)/2)*D_{i}/dxsq + hr3*D_{i+1})/dxsq = b
  !!
  !!
  !! 

  zero_out_linsys = .false.

  dx2 = 2.d0*dx
  dxsq = dx**2
  a3 = -alpha/3.d0

  DO i = 1,mx
      !! xc = xlower + (i-0.5)*dx
      hl = q(1,i-1)
      hc = q(1,i)
      hr = q(1,i+1)

      hl3 = ((hc+hl)/2.d0)**3
      hr3 = ((hc+hr)/2.d0)**3

      !!call bathy_complete(xc,b,dxzb,dx2zb)
      !!dxzb = bathy_deriv(xc)    !! Assume flat bathymetry for now

      dxzb  = aux(2,i)
      d2xzb = aux(3,i)
      dxh = (hr - hl)/dx2
      dxeta = dxh + dxzb

      bterm = alpha*(dxeta*dxzb + hc/2.d0*d2xzb)

      dm1 = a3*hl3/dxsq
      d0  = hc*(bterm + 1) - a3*(hr3 + hl3)/dxsq 
      dp1 = a3*hr3/dxsq

      IF (i .EQ. 1) THEN
          !! row 1
          T10   = dm1
          d(i)  = d0 + T10      !! + for Neumann; - for Dirichlet
          du(i) = dp1
      ELSE IF (i .EQ. mx) THEN
          !! row mx
          TNNp1   = dp1
          dl(i-1) = dm1
          d(i)    = d0 + TNNp1  !! + for Neumann; - for Dirichlet
      ELSE
          dl(i-1) = dm1
          d(i)    = d0
          du(i)   = dp1
      ENDIF    

  END DO

  !! Set up the right hand side
  CALL compute_rhs(meqn,mbc,mx,xlower,dx,q,maux,aux,b)

!! Zero out entries that do not pass test
  if (zero_out_linsys) then
      do i = 1,mx
          hc = q(1,i)

          dxzb  = aux(2,i)
          dxh = (hr - hl)/dx2
          !!dxh = (hc - hl)/dx
          dxeta = dxh + dxzb


          if (abs(dxeta) .gt. breaking .or. hc .le. dry_tolerance) then
              if (i .eq. 1) then
                  d(i) = 1
                  du(i) = 0.
              elseif (i .eq. mx) then
                  d(i) = 1.d0
                  dl(i-1) = 0.d0
              else
                  dl(i-1) = 0
                  d(i) = 1
                  du(i) = 0
              endif
                  b(i) = 0
          endif
      end do
  endif

  !! Solve tridiagonal system; solution stored in sgn_soln
  nrhs = 1
  ldb = mx
  CALL dgtsv(mx,nrhs,dl,d,du,b,ldb,info)
  !!call tri_solve(mx,dl,d,du,b)

  if (info .gt. 0) then
      write(6,*) 'info > 0; info = ', info
      stop
  endif

  do i = 1,mx
      sgn_soln(i) = b(i)
  end do

  do ibc = 1,mbc
      sgn_soln(1-ibc) = sgn_soln(ibc)
      sgn_soln(mx+ibc) = sgn_soln(mx-ibc+1)
  end do

  RETURN
END SUBROUTINE sgn

 
SUBROUTINE compute_rhs(meqn,mbc,mx,xlower,dx,q,maux,aux,b)  
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, maux
  DOUBLE PRECISION xlower, dx

  DOUBLE PRECISION  q(meqn,1-mbc:mx+mbc)
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc)
  DOUBLE PRECISION  b(mx)

  DOUBLE PRECISION :: grav, dry_tolerance, sea_level
  COMMON /common_swe/ grav, dry_tolerance, sea_level

  DOUBLE PRECISION :: alpha, breaking
  COMMON /common_dispersive/ alpha, breaking

  DOUBLE PRECISION  R1, R2
  DOUBLE PRECISION  w1(1-mbc:mx+mbc)
  DOUBLE PRECISION  w2(1-mbc:mx+mbc)

  DOUBLE PRECISION hc, hl, hr, dx2, dxw1, dxw2, dxh, dxzb
  DOUBLE PRECISION uc, ul, ur, divu, dxeta, d2xzb

  INTEGER i

  dx2 = 2*dx

  DO i = 2-mbc,mx+1

!!  double dxux = dx(u.x), dyuy = dy(u.y);
!!  c[] = - dxux*dyuy + dx(u.y)*dy(u.x) + sq(dxux + dyuy);
!!  d[] = sq(u.x[])*d2x(zb) + sq(u.y[])*d2y(zb) + 2.*u.x[]*u.y[]*d2xy(zb);
!!
!!        b.x[] = h[]*(G/alpha_d*dx(η) - 2.*R1(h,zb,c) + R2(h,zb,d));

    hl = q(1,i-1)
    ul = q(2,i-1)/hl

    hr = q(1,i+1)
    ur = q(2,i+1)/hr

    hc = q(1,i)
    uc = q(2,i)/hc

    divu = (ur - ul)/dx2

    d2xzb = aux(3,i)

    w1(i) =  divu**2       !! For 1d
    w2(i) = uc**2*d2xzb    !! For 1d
  END DO

  DO i = 1,mx
      hl = q(1,i-1)
      hc = q(1,i)
      hr = q(1,i+1)      

      dxh = (hr - hl)/dx2
      dxzb = aux(2,i)    !! Assume flat bathymetry for now
      dxeta = dxh + dxzb

      dxw1 = (w1(i+1) - w1(i-1))/dx2
      dxw2 = (w2(i+1) - w2(i-1))/dx2

      !! #define R1(h,zb,w) -h[]*(h[]/3.*dx(w) + w[]*(dx(h) + dx(zb)/2.))
      R1 = -hc*(hc/3.d0*dxw1 + w1(i)*(dxh + 0.5d0*dxzb))


      !! #define R2(h,zb,w) h[]/2.*dx(w) + w[]*(dx(zb) + dx(h))
      R2 = hc/2.d0*dxw2 + w2(i)*dxeta

!!    # Right hand side.  Note factor of h
!!        b.x[] = h[]*(G/alpha_d*dx(η) - 2.*R1(h,zb,c) + R2(h,zb,d));
      b(i) = hc*(grav/alpha*dxeta + (-2.d0*R1 + R2))
  END DO
!!  do i = mx-5, mx
!!    write(6,*) i, b(i)
!!  end do

END SUBROUTINE compute_rhs
