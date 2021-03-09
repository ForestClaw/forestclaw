SUBROUTINE src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
  IMPLICIT NONE

  INTEGER meqn, mbc, mx, maux
  DOUBLE PRECISION xlower, dx, t, dt

  DOUBLE PRECISION  q(meqn,1-mbc:mx+mbc)
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc)

  DOUBLE PRECISION :: grav, dry_tolerance, sea_level
  COMMON /common_swe/ grav, dry_tolerance, sea_level

  DOUBLE PRECISION :: alpha, breaking
  COMMON /common_dispersive/ alpha, breaking

  INTEGER i
  DOUBLE PRECISION hl, hc, hr, dxzb, dxh, dxeta, dx2
  !! DOUBLE PRECISION bathy

  DOUBLE PRECISION sgn_soln(1-mbc:mx+mbc), disp, xc, wall
  double precision wall_distance

  if (breaking .lt. 0) then
      return
  endif

  CALL sgn(meqn,mbc,mx,xlower,dx,q,maux,aux,sgn_soln)

  dx2 = 2.d0*dx

  DO i = 1,mx    !! Avoid problems at the wall?       
      !! xc = xlower + (i-0.5)*dx

      hr = q(1,i+1)
      hc = q(1,i)
      hl = q(1,i-1)


      dxzb = aux(2,i)            !! Gradient of the bathymetry 
      dxh = (hr - hl)/dx2
      dxeta = dxh + dxzb


      !!if (.true.) then
      if (abs(dxeta) .le. breaking .and. hc .ge. dry_tolerance) then
          disp = hc*(grav/alpha*dxeta - sgn_soln(i))
          q(2,i) = q(2,i) + dt*disp
      endif

  ENDDO


!!  sum(2) = 0
!!  do i = 1,mx
!!      sum(2) = sum(2) + q(2,i)*dx
!!  end do

!!  write(6,100) 'Before   : ', sum(1)
!!  write(6,100) 'After    : ', sum(2)
!!  if (sum(1) .ne. 0) then
!!    write(6,110) 'Diff (%) : ', 100*abs(sum(1)-sum(2))/abs(sum(1))
!!  else
!!    write(6,110) 'Diff (%) : ', abs(sum(1)-sum(2))
!!  endif    
!!  write(6,*)
!!
!!100 format(A12,E24.16)
!!110 format(A12,F8.2)


  RETURN
END SUBROUTINE src1
