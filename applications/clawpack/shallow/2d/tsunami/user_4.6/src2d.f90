SUBROUTINE src2d(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux,t,dt)
    IMPLICIT NONE

    INTEGER maxmx, maxmy, meqn, mbc, mx, my, maux
    DOUBLE PRECISION xlower, ylower, dx, dy, t, dt

    DOUBLE PRECISION  q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
    DOUBLE PRECISION  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    DOUBLE PRECISION :: alpha, breaking
    COMMON /common_dispersive/ alpha, breaking

    INTEGER i
    DOUBLE PRECISION hl, hc, hr, dxzb, dxh, dxeta, dx2
    DOUBLE PRECISION bathy

    DOUBLE PRECISION sgn_soln(1-mbc:mx+mbc), disp, xc, wall
    double precision wall_distance

    if (breaking .lt. 0) then
        return
    endif

    CALL sgn(meqn,mbc,mx,xlower,dx,q,maux,aux,sgn_soln)

    dx2 = 2.d0*dx

    DO i = 1,mx    !! Avoid problems at the wall?       
        !! xc = xlower + (i-0.5)*dx

        do j = 1,my
            hr = q(i+1,j,1)
            hc = q(i,j,1)
            hl = q(i-1,j,1)

            dxzb = aux(i,j,2)            !! Gradient of the bathymetry 
            dxh = (hr - hl)/dx2
            dxeta = dxh + dxzb

            !!if (.true.) then
            if (abs(dxeta) .le. breaking .and. hc .ge. dry_tolerance) then
              disp = hc*(grav/alpha*dxeta - sgn_soln(i))
              q(i,j,2) = q(i,j,2) + dt*disp
          endif
      end do

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
    END SUBROUTINE src2d
